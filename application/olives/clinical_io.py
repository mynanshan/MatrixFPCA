"""Read the OLIVES clinical workbook without an Excel-engine dependency."""
import csv
import re
import zipfile
import xml.etree.ElementTree as ET

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def read_xlsx(path, sheet=1):
    with zipfile.ZipFile(path) as z:
        ss = []
        if "xl/sharedStrings.xml" in z.namelist():
            ss = ["".join(e.itertext()) for e in ET.fromstring(
                z.read("xl/sharedStrings.xml")).findall("m:si", NS)]
        rows = []
        for _, r in ET.iterparse(z.open(f"xl/worksheets/sheet{sheet}.xml"), events=["end"]):
            if r.tag != "{" + NS["m"] + "}row":
                continue
            row = {}
            for c in r.findall("m:c", NS):
                v = c.find("m:v", NS)
                s = v.text if v is not None else ""
                if c.get("t") == "s":
                    s = ss[int(s)] if s else ""
                elif c.get("t") == "inlineStr":
                    s = "".join(c.find("m:is", NS).itertext())
                if s:
                    row[re.sub("[0-9]", "", c.get("r"))] = s
            rows.append(row)
            r.clear()
    return rows


def write_csv(path, records):
    if records:
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(records[0]))
            w.writeheader()
            w.writerows(records)

