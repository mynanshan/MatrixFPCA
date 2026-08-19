#!/usr/bin/env python3
"""Resumable UN Comtrade downloader for the MatrixFPCA trade application.

The API key is read from UN_COMTRADE_API_KEY. Each successful reporter/period
batch is written immediately, with a JSON sidecar used as a completion marker.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence


API_KEY_ENV = "UN_COMTRADE_API_KEY"
MAX_RECORDS = 100_000
DEFAULT_PERIOD_BATCH = 12
AVAILABILITY_CACHE_VERSION = 1


@dataclass(frozen=True)
class Country:
    region: str
    name: str
    iso3: str
    code: int


COUNTRIES: tuple[Country, ...] = (
    Country("North America", "United States", "USA", 842),
    Country("North America", "Canada", "CAN", 124),
    Country("North America", "Mexico", "MEX", 484),
    Country("Europe", "Germany", "DEU", 276),
    Country("Europe", "France", "FRA", 251),
    Country("Europe", "United Kingdom", "GBR", 826),
    Country("Europe", "Italy", "ITA", 380),
    Country("Europe", "Spain", "ESP", 724),
    Country("Europe", "Netherlands", "NLD", 528),
    Country("Europe", "Belgium", "BEL", 56),
    Country("Europe", "Switzerland", "CHE", 757),
    Country("Europe", "Sweden", "SWE", 752),
    Country("Europe", "Poland", "POL", 616),
    Country("East Asia", "China", "CHN", 156),
    Country("East Asia", "Japan", "JPN", 392),
    Country("East Asia", "South Korea", "KOR", 410),
    Country("East Asia", "Hong Kong SAR", "HKG", 344),
    Country("East Asia", "Singapore", "SGP", 702),
    Country("Other major economies", "India", "IND", 699),
    Country("Other major economies", "Brazil", "BRA", 76),
    Country("Other major economies", "Australia", "AUS", 36),
    Country("Other major economies", "T\u00fcrkiye", "TUR", 792),
    Country("Other major economies", "Saudi Arabia", "SAU", 682),
    Country("Other major economies", "South Africa", "ZAF", 710),
    Country("Other major economies", "Indonesia", "IDN", 360),
)

# HS chapter 77 is reserved and is not present in the requested specification.
HS2_CODES: tuple[str, ...] = tuple(
    f"{code:02d}" for code in range(1, 98) if code != 77
)


@dataclass(frozen=True)
class Query:
    reporter: Country
    periods: tuple[str, ...]

    @property
    def expected_cell_bound(self) -> int:
        return len(self.periods) * len(COUNTRIES) * len(HS2_CODES)

    @property
    def stem(self) -> str:
        return (
            f"reporter-{self.reporter.code:03d}_{self.reporter.iso3}_"
            f"periods-{self.periods[0]}-{self.periods[-1]}"
        )


class CallBudgetExhausted(RuntimeError):
    """Raised when the caller-requested API call budget has been used."""


class CallBudget:
    def __init__(self, maximum: int | None, min_interval: float = 0.0) -> None:
        self.maximum = maximum
        self.min_interval = min_interval
        self.used = 0
        self.last_started = 0.0

    def consume(self) -> None:
        if self.maximum is not None and self.used >= self.maximum:
            raise CallBudgetExhausted
        elapsed = time.monotonic() - self.last_started
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self.used += 1
        self.last_started = time.monotonic()


def month_range(start: str, end: str) -> tuple[str, ...]:
    """Return inclusive YYYYMM periods from YYYY-MM inputs."""

    def parse(value: str) -> tuple[int, int]:
        try:
            date = datetime.strptime(value, "%Y-%m")
        except ValueError as exc:
            raise ValueError(f"invalid month {value!r}; expected YYYY-MM") from exc
        return date.year, date.month

    start_year, start_month = parse(start)
    end_year, end_month = parse(end)
    if (start_year, start_month) > (end_year, end_month):
        raise ValueError("start month must not be after end month")

    result: list[str] = []
    year, month = start_year, start_month
    while (year, month) <= (end_year, end_month):
        result.append(f"{year:04d}{month:02d}")
        month += 1
        if month == 13:
            year += 1
            month = 1
    return tuple(result)


def batched(values: Sequence[str], size: int) -> Iterable[tuple[str, ...]]:
    if size < 1 or size > 12:
        raise ValueError("period batch size must be between 1 and 12")
    for offset in range(0, len(values), size):
        yield tuple(values[offset : offset + size])


def build_queries(
    start: str,
    end: str,
    reporters: Sequence[Country] = COUNTRIES,
    period_batch: int = DEFAULT_PERIOD_BATCH,
) -> list[Query]:
    periods = month_range(start, end)
    return [
        Query(reporter=reporter, periods=periods_batch)
        for periods_batch in batched(periods, period_batch)
        for reporter in reporters
    ]


def sdk_kwargs(query: Query) -> dict[str, Any]:
    """Translate the requested UI selections to comtradeapicall arguments."""

    return {
        "typeCode": "C",
        "freqCode": "M",
        "clCode": "HS",
        "period": ",".join(query.periods),
        "reporterCode": str(query.reporter.code),
        "cmdCode": ",".join(HS2_CODES),
        "flowCode": "X",
        "partnerCode": ",".join(str(country.code) for country in COUNTRIES),
        # Code 0/C00 means the total (all or unspecified) category. Supplying
        # these explicitly prevents accidental extended-detail rows.
        "partner2Code": "0",
        "customsCode": "C00",
        "motCode": "0",
        "maxRecords": MAX_RECORDS,
        "format_output": "JSON",  # the SDK only supports JSON -> DataFrame
        "aggregateBy": None,
        "breakdownMode": "classic",
        "countOnly": False,
        "includeDesc": True,
    }


def availability_kwargs(periods: Sequence[str]) -> dict[str, Any]:
    return {
        "typeCode": "C",
        "freqCode": "M",
        "clCode": "HS",
        "period": ",".join(periods),
        "reporterCode": ",".join(str(country.code) for country in COUNTRIES),
    }


def import_sdk() -> Any:
    try:
        import comtradeapicall
    except ImportError as exc:
        raise RuntimeError(
            "comtradeapicall is not installed; run: "
            "python -m pip install -r application/comtrade/requirements.txt"
        ) from exc
    return comtradeapicall


def frame_columns(frame: Any) -> set[str]:
    return {str(column) for column in frame.columns}


def validate_frame(frame: Any, query: Query) -> None:
    if frame is None:
        raise RuntimeError("the SDK returned no DataFrame")
    if not hasattr(frame, "columns") or not hasattr(frame, "to_csv"):
        raise RuntimeError("the SDK returned an unexpected object")
    if len(frame) == 0:
        raise RuntimeError("response contains no rows")
    required = {"period", "reporterCode", "partnerCode", "cmdCode", "primaryValue"}
    missing = required - frame_columns(frame)
    if missing:
        raise RuntimeError(f"response is missing columns: {', '.join(sorted(missing))}")
    if len(frame) >= MAX_RECORDS:
        raise RuntimeError(
            f"response reached maxRecords={MAX_RECORDS:,}; refusing a possibly truncated file"
        )

    def normalized(column: str, width: int = 0) -> set[str]:
        values = frame[column].dropna().astype(str).str.replace(r"\.0$", "", regex=True)
        if width:
            values = values.str.zfill(width)
        return set(values)

    checks = (
        ("period", normalized("period", 6), set(query.periods)),
        ("reporterCode", normalized("reporterCode"), {str(query.reporter.code)}),
        ("partnerCode", normalized("partnerCode"), {str(c.code) for c in COUNTRIES}),
        ("cmdCode", normalized("cmdCode", 2), set(HS2_CODES)),
    )
    for label, actual, allowed in checks:
        unexpected = actual - allowed
        if unexpected:
            sample = ", ".join(sorted(unexpected)[:5])
            raise RuntimeError(f"response has unexpected {label} value(s): {sample}")

    if "flowCode" in frame_columns(frame):
        unexpected_flows = normalized("flowCode") - {"X"}
        if unexpected_flows:
            raise RuntimeError(f"response has unexpected flows: {sorted(unexpected_flows)}")

    keys = ["period", "reporterCode", "partnerCode", "cmdCode"]
    duplicate_count = int(frame.duplicated(subset=keys).sum())
    if duplicate_count:
        raise RuntimeError(
            f"classic response has {duplicate_count:,} duplicate reporter-partner-HS2-month rows"
        )


def completion_paths(output_dir: Path, query: Query) -> tuple[Path, Path]:
    reporter_dir = output_dir / f"reporter={query.reporter.code:03d}"
    return reporter_dir / f"{query.stem}.csv", reporter_dir / f"{query.stem}.json"


def is_complete(
    csv_path: Path,
    metadata_path: Path,
    query: Query,
    available_periods: Sequence[str] | None = None,
) -> bool:
    if not metadata_path.is_file():
        return False
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    selection_parameters = metadata.get("selection_parameters", metadata.get("parameters"))
    if available_periods is not None:
        recorded_available = metadata.get("available_periods")
        if recorded_available is None:
            # Checkpoints made before availability preflighting are safe only
            # when every requested month was available.
            if list(available_periods) != list(query.periods):
                return False
        elif recorded_available != list(available_periods):
            return False
    common = (
        metadata.get("reporter_code") == query.reporter.code
        and metadata.get("periods") == list(query.periods)
        and selection_parameters == sdk_kwargs(query)
    )
    if metadata.get("status") == "unavailable":
        return common and metadata.get("row_count") == 0
    return (
        common
        and metadata.get("status") == "complete"
        and metadata.get("row_count", 0) > 0
        and csv_path.is_file()
        and csv_path.stat().st_size > 0
    )


def write_completed(
    frame: Any,
    output_dir: Path,
    query: Query,
    fetched_query: Query | None = None,
) -> Path:
    fetched_query = fetched_query or query
    csv_path, metadata_path = completion_paths(output_dir, query)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    csv_part = csv_path.with_suffix(".csv.part")
    metadata_part = metadata_path.with_suffix(".json.part")

    frame.to_csv(csv_part, index=False)
    metadata = {
        "status": "complete",
        "fetched_at_utc": datetime.now(timezone.utc).isoformat(),
        "reporter_code": query.reporter.code,
        "reporter_iso3": query.reporter.iso3,
        "periods": list(query.periods),
        "available_periods": list(fetched_query.periods),
        "unavailable_periods": [
            period for period in query.periods if period not in fetched_query.periods
        ],
        "row_count": int(len(frame)),
        "expected_cell_bound": query.expected_cell_bound,
        "selection_parameters": sdk_kwargs(query),
        "parameters": sdk_kwargs(fetched_query),
    }
    # Never put the subscription key in parameters or metadata.
    metadata_part.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    os.replace(csv_part, csv_path)
    os.replace(metadata_part, metadata_path)
    return csv_path


def write_unavailable(output_dir: Path, query: Query) -> Path:
    """Write a resolved checkpoint for a reporter-period batch with no monthly data."""

    _, metadata_path = completion_paths(output_dir, query)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    part = metadata_path.with_suffix(".json.part")
    metadata = {
        "status": "unavailable",
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "reporter_code": query.reporter.code,
        "reporter_iso3": query.reporter.iso3,
        "periods": list(query.periods),
        "available_periods": [],
        "unavailable_periods": list(query.periods),
        "row_count": 0,
        "expected_cell_bound": query.expected_cell_bound,
        "selection_parameters": sdk_kwargs(query),
        "reason": "No monthly HS dataset is listed by UN Comtrade data availability.",
    }
    part.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(part, metadata_path)
    return metadata_path


def fetch_with_retries(
    fetch: Callable[..., Any],
    api_key: str,
    query: Query,
    budget: CallBudget,
    retries: int,
) -> Any:
    for attempt in range(retries + 1):
        budget.consume()
        try:
            frame = fetch(api_key, **sdk_kwargs(query))
            validate_frame(frame, query)
            return frame
        except CallBudgetExhausted:
            raise
        except Exception as exc:
            if attempt == retries:
                raise RuntimeError(
                    f"failed {query.reporter.iso3} {query.periods[0]}-{query.periods[-1]} "
                    f"after {retries + 1} attempts"
                ) from exc
            delay = min(60.0, 2.0**attempt) + random.uniform(0.0, 0.5)
            print(f"  retry {attempt + 1}/{retries} after {delay:.1f}s: {exc}", file=sys.stderr)
            time.sleep(delay)
    raise AssertionError("unreachable")


def load_availability_cache(
    output_dir: Path, periods: Sequence[str]
) -> set[tuple[int, str]] | None:
    path = output_dir / "availability.json"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if (
        payload.get("version") != AVAILABILITY_CACHE_VERSION
        or payload.get("periods") != list(periods)
        or payload.get("reporter_codes") != [country.code for country in COUNTRIES]
    ):
        return None
    try:
        return {
            (int(item["reporter_code"]), str(item["period"]))
            for item in payload["available_reporter_periods"]
        }
    except (KeyError, TypeError, ValueError):
        return None


def fetch_availability(
    sdk: Any,
    output_dir: Path,
    periods: Sequence[str],
    period_batch: int,
    budget: CallBudget,
    retries: int,
    overwrite: bool,
) -> set[tuple[int, str]]:
    """Fetch and cache monthly dataset availability using the public endpoint."""

    if not overwrite:
        cached = load_availability_cache(output_dir, periods)
        if cached is not None:
            print(f"Loaded cached data availability ({len(cached):,} reporter-months).")
            return cached

    available: set[tuple[int, str]] = set()
    batches = list(batched(periods, period_batch))
    for index, periods_batch in enumerate(batches, start=1):
        label = f"availability {periods_batch[0]}-{periods_batch[-1]}"
        for attempt in range(retries + 1):
            budget.consume()
            try:
                frame = sdk._getFinalDataAvailability(**availability_kwargs(periods_batch))
                if frame is None or not hasattr(frame, "columns"):
                    raise RuntimeError("availability API returned no DataFrame")
                required = {"period", "reporterCode"}
                missing = required - frame_columns(frame)
                if missing:
                    raise RuntimeError(
                        f"availability response is missing columns: {', '.join(sorted(missing))}"
                    )
                if len(frame) == 0:
                    raise RuntimeError("availability response contains no rows for all 25 reporters")
                for reporter_code, period in zip(frame["reporterCode"], frame["period"]):
                    pair = (int(reporter_code), str(period).replace(".0", "").zfill(6))
                    if pair[0] in {country.code for country in COUNTRIES} and pair[1] in periods_batch:
                        available.add(pair)
                print(
                    f"[{index}/{len(batches)}] {label}: "
                    f"{sum(1 for pair in available if pair[1] in periods_batch)} reporter-months"
                )
                break
            except Exception as exc:
                if attempt == retries:
                    raise RuntimeError(f"failed {label} after {retries + 1} attempts") from exc
                delay = min(60.0, 2.0**attempt) + random.uniform(0.0, 0.5)
                print(
                    f"  availability retry {attempt + 1}/{retries} after {delay:.1f}s: {exc}",
                    file=sys.stderr,
                )
                time.sleep(delay)

    payload = {
        "version": AVAILABILITY_CACHE_VERSION,
        "fetched_at_utc": datetime.now(timezone.utc).isoformat(),
        "periods": list(periods),
        "reporter_codes": [country.code for country in COUNTRIES],
        "available_reporter_periods": [
            {"reporter_code": reporter_code, "period": period}
            for reporter_code, period in sorted(available)
        ],
    }
    path = output_dir / "availability.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(".json.part")
    part.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)
    return available


def write_manifest(output_dir: Path, queries: Sequence[Query], args: argparse.Namespace) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "manifest.json"
    payload = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "start": args.start,
        "end": args.end,
        "period_batch": args.period_batch,
        "country_count": len(COUNTRIES),
        "hs2_count": len(HS2_CODES),
        "query_count": len(queries),
        "checkpoint_csvs": [
            str(completion_paths(output_dir, query)[0].relative_to(output_dir))
            for query in queries
        ],
        "countries": [asdict(country) for country in COUNTRIES],
        "hs2_codes": list(HS2_CODES),
        "api_key_environment_variable": API_KEY_ENV,
    }
    part = path.with_suffix(".json.part")
    part.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(part, path)


def show_plan(queries: Sequence[Query], start: str, end: str) -> None:
    periods = month_range(start, end)
    dense_cells = len(COUNTRIES) ** 2 * len(HS2_CODES) * len(periods)
    print(f"Countries (reporters and partners): {len(COUNTRIES)}")
    print(f"HS2 codes: {len(HS2_CODES)} (01-76 and 78-97; 77 excluded)")
    print(f"Periods: {len(periods)} ({periods[0]}-{periods[-1]})")
    print(f"Dense target cells: {dense_cells:,}")
    availability_calls = (len(periods) + DEFAULT_PERIOD_BATCH - 1) // DEFAULT_PERIOD_BATCH
    print(f"Trade-data queries (maximum): {len(queries)}")
    print(f"Data-availability queries (maximum, then cached): {availability_calls}")
    print(f"Largest requested query: {max(q.expected_cell_bound for q in queries):,} cells")
    print(f"Per-call record cap: {MAX_RECORDS:,}")


def run_download(args: argparse.Namespace, queries: Sequence[Query]) -> int:
    api_key = os.environ.get(API_KEY_ENV, "").strip()
    if not api_key:
        print(
            f"Missing {API_KEY_ENV}. Export the rotated key before running the downloader.",
            file=sys.stderr,
        )
        return 2

    output_dir = Path(args.output_dir).expanduser().resolve()
    write_manifest(output_dir, queries, args)
    sdk = import_sdk()
    budget = CallBudget(args.max_calls, min_interval=args.min_interval)
    try:
        availability = fetch_availability(
            sdk,
            output_dir,
            month_range(args.start, args.end),
            period_batch=args.period_batch,
            budget=budget,
            retries=args.retries,
            overwrite=args.overwrite,
        )
    except CallBudgetExhausted:
        print(
            f"Stopped cleanly after {budget.used} API calls while checking availability; "
            "rerun the same command to resume.",
            file=sys.stderr,
        )
        return 3

    completed = skipped = unavailable = 0

    for index, query in enumerate(queries, start=1):
        csv_path, metadata_path = completion_paths(output_dir, query)
        label = f"[{index}/{len(queries)}] {query.reporter.iso3} {query.periods[0]}-{query.periods[-1]}"
        available_periods = tuple(
            period
            for period in query.periods
            if (query.reporter.code, period) in availability
        )
        if (
            is_complete(
                csv_path,
                metadata_path,
                query,
                available_periods=available_periods,
            )
            and not args.overwrite
        ):
            skipped += 1
            print(f"{label}: already complete")
            continue
        if not available_periods:
            path = write_unavailable(output_dir, query)
            unavailable += len(query.periods)
            print(f"{label}: no monthly dataset; recorded gap in {path}")
            continue
        fetched_query = Query(reporter=query.reporter, periods=available_periods)
        if len(available_periods) != len(query.periods):
            missing_count = len(query.periods) - len(available_periods)
            unavailable += missing_count
            print(f"{label}: downloading {len(available_periods)} available months ({missing_count} unavailable)")
        else:
            print(f"{label}: downloading")
        try:
            frame = fetch_with_retries(
                sdk.getFinalData,
                api_key,
                fetched_query,
                budget,
                retries=args.retries,
            )
        except CallBudgetExhausted:
            print(
                f"Stopped cleanly after {budget.used} API calls; rerun the same command to resume.",
                file=sys.stderr,
            )
            return 3
        path = write_completed(frame, output_dir, query, fetched_query=fetched_query)
        completed += 1
        print(f"  wrote {len(frame):,} rows to {path}")

    print(
        f"Finished: {completed} downloaded, {skipped} already complete, "
        f"{unavailable} unavailable reporter-months, {budget.used} API calls this run."
    )
    return 0


def run_combine(args: argparse.Namespace) -> int:
    """Stream completed checkpoints into one CSV without loading them all."""

    output_dir = Path(args.output_dir).expanduser().resolve()
    destination = Path(args.combined_file).expanduser()
    if not destination.is_absolute():
        destination = output_dir / destination
    destination = destination.resolve()
    if destination.exists() and not args.overwrite:
        print(f"Refusing to replace {destination}; pass --overwrite to rebuild it.", file=sys.stderr)
        return 2

    manifest_path = output_dir / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        checkpoint_names = manifest["checkpoint_csvs"]
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as exc:
        print(f"Cannot read a valid checkpoint list from {manifest_path}: {exc}", file=sys.stderr)
        return 2
    sources = [output_dir / name for name in checkpoint_names]
    complete_sources: list[Path] = []
    resolved_count = 0
    availability_gaps: list[dict[str, Any]] = []
    for source in sources:
        metadata_path = source.with_suffix(".json")
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        status = metadata.get("status")
        if status == "complete" and metadata.get("row_count", 0) > 0 and source.is_file():
            complete_sources.append(source)
            resolved_count += 1
        elif status == "unavailable":
            resolved_count += 1
        for period in metadata.get("unavailable_periods", []):
            availability_gaps.append(
                {
                    "reporterCode": metadata.get("reporter_code"),
                    "reporterISO": metadata.get("reporter_iso3"),
                    "period": period,
                    "reason": metadata.get(
                        "reason", "No monthly HS dataset is listed by UN Comtrade data availability."
                    ),
                }
            )
    if not complete_sources:
        print(f"No completed checkpoint CSVs found below {output_dir}.", file=sys.stderr)
        return 2
    if resolved_count != len(sources) and not args.allow_partial:
        print(
            f"Only {resolved_count} of {len(sources)} manifest checkpoints are resolved; "
            "finish the download or pass --allow-partial intentionally.",
            file=sys.stderr,
        )
        return 2

    destination.parent.mkdir(parents=True, exist_ok=True)
    part = destination.with_suffix(destination.suffix + ".part")
    header: list[str] = []
    seen_columns: set[str] = set()
    for source in complete_sources:
        with source.open("r", encoding="utf-8", newline="") as handle:
            source_header = next(csv.reader(handle), None)
        if not source_header:
            raise RuntimeError(f"empty checkpoint: {source}")
        for column in source_header:
            if column not in seen_columns:
                header.append(column)
                seen_columns.add(column)
    if "primaryValue" not in seen_columns:
        raise RuntimeError("checkpoint CSVs do not contain primaryValue")

    row_count = 0
    try:
        with part.open("w", encoding="utf-8", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=header)
            writer.writeheader()
            for source in complete_sources:
                with source.open("r", encoding="utf-8", newline="") as handle:
                    reader = csv.DictReader(handle)
                    for row in reader:
                        if None in row:
                            raise RuntimeError(f"malformed CSV row in checkpoint: {source}")
                        writer.writerow(row)
                        row_count += 1
        os.replace(part, destination)
    except Exception:
        try:
            part.unlink()
        except FileNotFoundError:
            pass
        raise
    gaps_path = destination.with_name(f"{destination.stem}_availability_gaps.csv")
    gaps_part = gaps_path.with_suffix(gaps_path.suffix + ".part")
    with gaps_part.open("w", encoding="utf-8", newline="") as handle:
        gap_writer = csv.DictWriter(
            handle, fieldnames=("reporterCode", "reporterISO", "period", "reason")
        )
        gap_writer.writeheader()
        gap_writer.writerows(availability_gaps)
    os.replace(gaps_part, gaps_path)
    print(
        f"Combined {len(complete_sources)} completed checkpoints and {row_count:,} rows into "
        f"{destination}"
    )
    print(f"Recorded {len(availability_gaps):,} unavailable reporter-months in {gaps_path}")
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "command",
        choices=("plan", "pilot", "download", "combine"),
        help=(
            "plan without network access, run a one-month pilot, download the requested range, "
            "or combine completed checkpoints"
        ),
    )
    result.add_argument("--start", default="2012-01", help="first month, YYYY-MM")
    result.add_argument("--end", default="2024-12", help="last month, YYYY-MM")
    result.add_argument(
        "--pilot-period", default="2024-01", help="one month used by the pilot command, YYYY-MM"
    )
    result.add_argument(
        "--output-dir", default="data/comtrade", help="checkpoint/output directory"
    )
    result.add_argument(
        "--combined-file",
        default="comtrade_hs2_exports.csv",
        help="combined CSV path (relative paths are below --output-dir)",
    )
    result.add_argument(
        "--period-batch", type=int, default=DEFAULT_PERIOD_BATCH, help="months per API query (1-12)"
    )
    result.add_argument(
        "--max-calls", type=int, default=None, help="stop cleanly after this many calls, including retries"
    )
    result.add_argument("--retries", type=int, default=4, help="retries after the first attempt")
    result.add_argument(
        "--min-interval", type=float, default=1.05, help="minimum seconds between attempts"
    )
    result.add_argument(
        "--overwrite", action="store_true", help="replace completed checkpoint files"
    )
    result.add_argument(
        "--allow-partial",
        action="store_true",
        help="allow combine to use fewer completed checkpoints than the manifest requires",
    )
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.max_calls is not None and args.max_calls < 1:
        raise SystemExit("--max-calls must be positive")
    if args.retries < 0:
        raise SystemExit("--retries must be nonnegative")
    if args.min_interval < 0:
        raise SystemExit("--min-interval must be nonnegative")
    if args.command == "combine":
        return run_combine(args)

    try:
        if args.command == "pilot":
            pilot = month_range(args.pilot_period, args.pilot_period)[0]
            args.start = f"{pilot[:4]}-{pilot[4:]}"
            args.end = args.start
            args.period_batch = 1
        queries = build_queries(args.start, args.end, period_batch=args.period_batch)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    show_plan(queries, args.start, args.end)
    if args.command == "plan":
        return 0
    return run_download(args, queries)


if __name__ == "__main__":
    raise SystemExit(main())
