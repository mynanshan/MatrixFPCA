# UN Comtrade monthly bilateral HS2 exports

This downloader implements the requested 25-reporter by 25-partner monthly
export panel for all HS2 chapters except reserved chapter 77, from January 2012
through December 2024. It uses UN Comtrade's final-data API through the official
`comtradeapicall` Python package.

## Security first

The keys visible in the supplied screenshot must be treated as compromised.
Regenerate both keys in the UN Comtrade developer portal. Do not paste the new
key into this script, a notebook, or Git. Provide it only through the
environment variable used below. `.env` files are ignored by this project, but
the downloader deliberately does not load them automatically.

## API mapping

| UI selection | API/package argument |
|---|---|
| Goods | `typeCode="C"` |
| Monthly | `freqCode="M"` |
| HS, As Reported | `clCode="HS"` |
| Exports | `flowCode="X"` |
| Second partner: All | `partner2Code="0"` |
| Customs: All | `customsCode="C00"` |
| Mode of transport: All | `motCode="0"` |
| Classic | `breakdownMode="classic"` |
| Aggregation: None | `aggregateBy=None` |

The package fetches JSON into a pandas DataFrame; the script then writes CSV,
including the API's `primaryValue` field. Passing `format_output="CSV"` to the
package is not used because its DataFrame function currently supports JSON
only.

The script uses current Comtrade trade-area codes rather than blindly copying
standard M49 codes. In particular, France is `251` (France and Monaco),
Switzerland is `757` (Switzerland and Liechtenstein), and India is `699`. The
codes `356` and `756` from the initial table are not valid current Comtrade
reporter codes: India `356` expired in 1974 and Switzerland `756` is partner-only.

## Setup

From the project root:

```bash
python3 -m venv .venv-comtrade
source .venv-comtrade/bin/activate
python -m pip install -r application/comtrade/requirements.txt
read -rsp 'UN Comtrade API key: ' UN_COMTRADE_API_KEY
export UN_COMTRADE_API_KEY
echo
```

First inspect the no-network request plan:

```bash
python application/comtrade/download_comtrade.py plan
```

Then run the recommended one-month sanity check. It downloads one checkpoint
for each of the 25 reporters, with the same 25 partners and all 96 HS2 codes:

```bash
python application/comtrade/download_comtrade.py pilot --pilot-period 2024-01 \
  --output-dir data/comtrade-pilot
```

Inspect a few pilot CSV files and confirm these fields:

- `flowCode` is `X`;
- `reporterCode` is the reporter named in the filename;
- `partnerCode` contains only the 25 requested areas;
- `cmdCode` contains only the requested two-digit chapters;
- `primaryValue` is present;
- the key `(period, reporterCode, partnerCode, cmdCode)` is unique.

The script checks all of these automatically except that not every requested
combination must have a row: Comtrade generally omits combinations with no
reported trade. Absence of a row should not automatically be interpreted as a
zero when the reporter-month dataset itself is unavailable.

Once the pilot looks right, run the complete extraction:

```bash
python application/comtrade/download_comtrade.py download \
  --start 2012-01 --end 2024-12 --output-dir data/comtrade
```

After all checkpoints are complete, stream them into one CSV without loading
the full result into memory:

```bash
python application/comtrade/download_comtrade.py combine \
  --output-dir data/comtrade
```

This creates `data/comtrade/comtrade_hs2_exports.csv`. Rebuilding an existing
combined file requires `--overwrite`. Combining is refused if any checkpoint
listed in the current manifest is missing; `--allow-partial` opts into an
incomplete combined file explicitly.

The full plan makes at most 325 trade-data calls: 25 reporters for each of 13
batches, where one batch contains at most 12 months. Before downloading, 13
public data-availability calls determine which reporter-month datasets actually
exist; this result is cached in `availability.json`. Unavailable months are
skipped and recorded rather than retried as errors. Each trade request has at
most 28,800 reporter-partner-HS2-month cells, below the free plan's
100,000-record limit.
The registered free plan currently permits 500 calls per day at one call per
second, so the default 1.05-second spacing is intentional. Retries also count
against the quota.

Every completed call creates a CSV and a JSON sidecar below a reporter-specific
directory. The sidecar is the completion marker and records the exact nonsecret
parameters, periods, row count, and retrieval time. Repeating the same command
skips completed files and resumes unfinished work. Partial `.part` files are
never considered complete.

Some reporters do not supply monthly data for every year even when annual data
exists. For example, UN Comtrade lists no monthly HS dataset for the Republic of
Korea in 2012. These are genuine source-data gaps, not zero trade. After
combining, `comtrade_hs2_exports_availability_gaps.csv` lists every unavailable
reporter-month. Treat those reporter-month matrix rows as missing (`NA`), not as
zeros.

To reserve some of the daily quota, use for example `--max-calls 450`. The
process exits cleanly when that budget is reached; rerun the same command later
to continue. Use `--overwrite` only when intentionally refreshing the latest
Comtrade revisions.

## Verification

The specification and request-plan tests require no API key or third-party
packages:

```bash
python -m unittest application/comtrade/test_download_comtrade.py
```

UN Comtrade serves the latest available version of a dataset and may revise it
later. Keep the JSON retrieval timestamps if reproducibility matters, and use
`--overwrite` for a deliberate refresh rather than silently mixing vintages.
