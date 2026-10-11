"""o5c lane B: src-generic-12 (Zenodo 22150094, 'example sEEG ... with clinical confirmed one seizure') -> iEEG-BIDS.
One EDF+ file (sEEG-HFOs-8.edf, 103 signals incl. EDF Annotations, 256 Hz, 14095 s) copied byte-identical when the
identification fields are already anonymised (they are placeholders 'patient_2_3_17 X 01-JAN-1999 name' / 01-JAN-2019);
otherwise the header is scrubbed (data records untouched). Annotations -> events.tsv. The post-implantation MRI (DICOM,
T1/T2 with electrodes, not defaced) is NOT included (privacy rule: T1w/CT must be defaced or withheld).
Usage: python o5c_B_conv_generic12.py src-generic-12"""
import json, os, re, sys, urllib.request
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o5c_B_common import *  # noqa

CID = sys.argv[1] if len(sys.argv) > 1 else "src-generic-12"
SRC = src_dir(CID); OUT = fresh_out(CID); REP = rep_dir(CID)
EDF = os.path.join(SRC, "files", "sEEG-HFOs-8.edf")
ZEN = "22150094"; DOI = f"10.5281/zenodo.{ZEN}"; URL = f"https://zenodo.org/records/{ZEN}"
CONV = "harness/o5c_B_conv_generic12.py (ieeg-nemar-20261005 campaign, lane o5c-B)"
rep = {"cand_id": CID, "converter": CONV, "source_edf": EDF}

# repository record (metadata only) -> sourcedata
try:
    rec = json.load(urllib.request.urlopen(f"https://zenodo.org/api/records/{ZEN}", timeout=60))
except Exception as e:
    rec = {"error": repr(e)}
sd = os.path.join(OUT, "sourcedata", f"zenodo-{ZEN}"); os.makedirs(sd, exist_ok=True)
wjson(os.path.join(sd, f"zenodo_record_{ZEN}.json"), rec)
meta = rec.get("metadata", {})
authors = [c["name"] for c in meta.get("creators", [])] or ["AI-is-out-there"]
lic = (meta.get("license") or {}).get("id", "cc-by-4.0").upper().replace("CC-BY-4.0", "CC-BY-4.0")
rep["licence"] = {"value": lic, "read_from": "zenodo api metadata.license.id"}
rep["record_description"] = re.sub("<[^>]+>", " ", meta.get("description", ""))

# EDF header review
hdr = edf_header(EDF)
rep["edf_header"] = {k: hdr[k] for k in ("patient", "recording", "startdate", "starttime", "n_records", "rec_dur", "ns", "reserved", "file_bytes", "expected_bytes")}
assert hdr["file_bytes"] == hdr["expected_bytes"], "EDF size mismatch"
ann, rec_on = edf_annotations(EDF, hdr, encoding="cp1251")  # texts are Russian (Windows-1251): "где тут начало?", "клиника", "приступ + БТКП"
rep["n_annotations"] = len(ann)
rep["annotation_texts_unique"] = sorted({a[2] for a in ann})[:300]
rep["annotation_flags"] = flag([a[2] for a in ann] + [hdr["patient"], hdr["recording"]])
rep["annotation_encoding"] = "cp1251 (Russian)"
# identification fields: placeholders? (no real name/date: 'X' and 01-JAN-1999 / 01-JAN-2019 style)
pat_ok = bool(re.fullmatch(r"patient_[\d_]+ X 01-JAN-\d{4} name", hdr["patient"])) and hdr["recording"].startswith("Startdate 01-JAN-")
rep["identification_fields_placeholders"] = pat_ok

sub = "sub-01"; task = "seizure"
D = os.path.join(OUT, sub, "ieeg"); os.makedirs(D, exist_ok=True)
stem = f"{sub}_task-{task}"
dst = os.path.join(D, stem + "_ieeg.edf")
if pat_ok and not rep["annotation_flags"]:
    cp = copy_verbatim(EDF, dst); rep["edf_copy"] = cp; rep["edf_mode"] = "byte-identical copy"
    assert cp["identical"]
else:
    cp = edf_scrub_header(EDF, dst, patient="X X X X"); rep["edf_copy"] = cp; rep["edf_mode"] = "header scrubbed (patient/recording/startdate day), data records identical"

# channels
labels = hdr["labels"]
keep = [i for i, l in enumerate(labels) if l.upper() != "EDF ANNOTATIONS"]
names = [labels[i] for i in keep]
def ctype(l):
    if l.startswith("EEG "):
        return "SEEG"
    if l.startswith("MKR"):
        return "MISC"
    return "MISC"
types = [ctype(l) for l in names]
units = [hdr["units"][i] or "n/a" for i in keep]
units = ["uV" if u in ("uV", "µV") else u for u in units]
def cut(pf):
    m1 = re.search(r"HP:([\d.]+)Hz", pf); m2 = re.search(r"LP:([\d.]+)Hz", pf)
    return (float(m1.group(1)) if m1 else "n/a", float(m2.group(1)) if m2 else "n/a")
lows = [cut(hdr["prefilter"][i])[0] for i in keep]; highs = [cut(hdr["prefilter"][i])[1] for i in keep]
def group(l):
    if not l.startswith("EEG "):
        return "n/a"
    return re.sub(r"\d+$", "", l[4:])
extra = [{"group": group(l), "description": ("sEEG contact (EDF label kept)" if t == "SEEG" else "event-marker input channel of the EEG system (EDF label kept)")}
         for l, t in zip(names, types)]
sf = hdr["sfreq"][keep[0]]
write_channels(os.path.join(D, stem + "_channels.tsv"), names, types, units, sf, extra=extra, low=lows, high=highs)
groups = sorted({group(l) for l in names if group(l) != "n/a"}, key=lambda g: names.index(next(n for n in names if group(n) == g)))
rep["channels"] = {"n": len(names), "seeg": types.count("SEEG"), "misc": types.count("MISC"), "electrodes": groups, "sfreq": sf}

# events from annotations
EN = {"где тут начало?": "where is the onset here? (reviewer's note)", "клиника": "clinical signs (semiology) begin", "приступ + БТКП": "seizure + BTCS (bilateral tonic-clonic seizure)"}
ev = [[round(a[0], 6), (round(a[1], 6) if a[1] is not None else "n/a"), int(round(a[0] * sf)), a[2], EN.get(a[2].strip(), "n/a")] for a in ann]
if ev:
    wtsv(os.path.join(D, stem + "_events.tsv"), ["onset", "duration", "sample", "trial_type", "trial_type_english"], ev)
    wjson(os.path.join(D, stem + "_events.json"), {"onset": {"Description": "EDF+ annotation onset, seconds from the start of the file", "Units": "s"},
          "duration": {"Description": "EDF+ annotation duration (n/a when the annotation has none)", "Units": "s"},
          "sample": {"Description": "onset x 256 Hz, rounded"},
          "trial_type": {"Description": "EDF+ annotation text, verbatim (Russian, stored in the file as Windows-1251 and decoded here to UTF-8; the release does not document the annotation vocabulary)"},
          "trial_type_english": {"Description": "English rendering of the annotation added by the converter (not part of the release); n/a if not translated"}})

wjson(os.path.join(D, stem + "_ieeg.json"), ieeg_sidecar(
    "seizure", sf, hdr["duration_s"], types, "n/a", "n/a",
    TaskDescription="Clinical stereo-EEG monitoring segment containing one clinically confirmed seizure (Zenodo record: \"seizure 17:27:11, clinicaly confirmed 17:27:14\"; the EDF clock fields are anonymised placeholders, so the seizure time is not recoverable from the header; see events.tsv for the annotations present in the file).",
    Instructions="n/a", InstitutionName="n/a", ElectrodeManufacturer="n/a", HardwareFilters={"EDF prefilter field, sEEG channels": {"HighpassHz": 0.15, "LowpassHz": 67}, "EDF prefilter field, marker channels": {"HighpassHz": 0, "LowpassHz": 67}},
    ElectricalStimulation=False, RecordingType="continuous"))
# electrodes (no coordinates released outside the withheld MRI)
wtsv(os.path.join(D, f"{sub}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size", "group"],
     [[n, "n/a", "n/a", "n/a", "n/a", group(n)] for n, t in zip(names, types) if t == "SEEG"])
wjson(os.path.join(D, f"{sub}_space-Other_coordsystem.json"), {"iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "n/a",
      "iEEGCoordinateSystemDescription": "No electrode coordinates are released as numbers. The Zenodo record contains a post-implantation MRI (DICOM T1/T2 with the electrodes visible) that is not defaced; it is withheld from this BIDS dataset for privacy and can be obtained from the Zenodo record (MRI-with-electrodes.zip). x, y, z = n/a."})
wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "age", "sex", "group"], [[sub, "n/a", "n/a", "epilepsy"]])
wjson(os.path.join(OUT, "participants.json"), {"age": {"Description": "not released"}, "sex": {"Description": "not released"},
                                               "group": {"Description": "clinical population stated by the record title (sEEG monitoring with one clinically confirmed seizure)"}})
name = "Zenodo-22150094: example stereo-EEG (sEEG) recording with one clinically confirmed seizure"
wjson(os.path.join(OUT, "dataset_description.json"), dataset_description(
    name, authors, "CC-BY-4.0", DOI, URL, [f"https://doi.org/{DOI}"],
    f"Cite the Zenodo record doi:{DOI} ({authors[0]}, version {meta.get('version', 'n/a')}, {meta.get('publication_date', 'n/a')}).",
    ["n/a: the Zenodo record does not state an ethics approval; the recording was shared by its depositor under CC-BY-4.0, see the record"],
    [{"Name": "o5c_B_conv_generic12.py", "Description": "Byte-identical copy of the released EDF+ file into iEEG-BIDS with channels/events/electrodes sidecars derived from its header and annotations; no signal processing. Code in code/."}],
    keywords=["stereo-EEG", "SEEG", "epilepsy", "seizure", "HFO"], source_version=f"Zenodo version {meta.get('version', 'n/a')} ({meta.get('publication_date', 'n/a')})"))
ann_lines = "\n".join(f"  - {a[0]:.3f} s: {a[2]}" for a in ann[:40]) or "  (none)"
open(os.path.join(OUT, "README"), "w").write(f"""# Example stereo-EEG recording with one clinically confirmed seizure (Zenodo 22150094)

Stereo-EEG (sEEG) recording of one patient during clinical intracranial monitoring, released on Zenodo by
"{authors[0]}" as "{meta.get('title', 'example sEEG ... with clinical confirmed one seizure')}" (doi:{DOI}, CC-BY-4.0,
version {meta.get('version', 'n/a')}). The record's whole description reads: "'sEEG-HFOs-8.edf': seizure 17:27:11,
clinicaly confirmed 17:27:14 / MRI T1 and T2". No article, protocol, ethics statement or demographic information is
given in the record; this README does not add any.

## Content
- `sub-01/ieeg/sub-01_task-seizure_ieeg.edf`: the released EDF+ file, {rep['edf_mode']} (sha256
  {cp.get('sha256', cp.get('sha256_dst'))}). {hdr['ns'] - 1} signals ({types.count('SEEG')} sEEG contacts labelled `EEG <electrode><contact>`,
  {types.count('MISC')} marker inputs `MKR1+`/`MKR2+`) at {sf:g} Hz, {hdr['duration_s'] / 3600:.2f} h ({hdr['n_records']} one-second records), physical range
  +-3200 uV, 16-bit. The EDF identification fields are the depositor's placeholders (`{hdr['patient']}`, `{hdr['recording']}`),
  so the recording date/time in the header is not the real one.
- Electrodes ({len(groups)}): {', '.join(groups)}. Names ending in `'` follow the usual sEEG convention for the
  contralateral hemisphere, but the record does not state which hemisphere is which, so no hemisphere is asserted here.
- `events.tsv`: the {len(ann)} EDF+ annotations of the file, verbatim (Russian, Windows-1251 in the file, UTF-8 here) with an
  English rendering added by the converter. They mark the clinical onset and a seizure with bilateral tonic-clonic
  evolution ("приступ + БТКП").
{ann_lines}
- Not included: `MRI-with-electrodes.zip` (DICOM T1 and T2 with implanted electrodes, {rep.get('dicom_files', 432)} slices). It is not
  defaced, so it is withheld from this redistribution; it is available from the Zenodo record. No electrode coordinates
  are therefore released (`electrodes.tsv` has n/a coordinates).
- `sourcedata/zenodo-{ZEN}/`: the Zenodo record metadata (JSON).

## Provenance
Converted by the ieeg-nemar-20261005 campaign (lane o5c-B) with `code/o5c_B_conv_generic12.py`; the EDF data records are
unchanged. `code/conversion_report.json` lists the header fields, the annotation texts and the checksums.
""")
copy_code(CID, ["o5c_B_conv_generic12.py", "o5c_B_common.py"])
bidsignore(OUT, ["sourcedata/", "code/"])

# round trip: MNE read of the BIDS file vs digital values from the source bytes
import mne
raw = mne.io.read_raw_edf(dst, preload=False, verbose="error", encoding="latin1")
chk = {}
for ci in (0, 10, len(keep) - 2):
    i = keep[ci]
    dig = edf_signal_check(EDF, hdr, i, n_rec=20)
    scale = (hdr["pmax"][i] - hdr["pmin"][i]) / (hdr["dmax"][i] - hdr["dmin"][i])
    phys = (dig - hdr["dmin"][i]) * scale + hdr["pmin"][i]
    x = raw.get_data(picks=[raw.ch_names.index(names[ci])], start=0, stop=len(dig))[0] * 1e6
    chk[names[ci]] = {"max_abs_diff_uV": float(np.max(np.abs(x - phys))), "n": int(len(dig))}
rep["roundtrip"] = {"mne_n_times": raw.n_times, "expected_n_times": int(hdr["n_records"] * hdr["nsamp"][keep[0]]), "mne_n_ch": len(raw.ch_names), "channels": chk,
                    "ok": raw.n_times == hdr["n_records"] * hdr["nsamp"][keep[0]] and all(v["max_abs_diff_uV"] < 1e-6 for v in chk.values())}
rep["hours"] = hdr["duration_s"] / 3600; rep["subjects"] = 1; rep["n_recordings"] = 1
report(CID, rep)
print(json.dumps({k: rep[k] for k in ("edf_mode", "n_annotations", "annotation_flags", "roundtrip", "hours")}, indent=1, default=str))
print("CONVERT_OK" if rep["roundtrip"]["ok"] else "CONVERT_ROUNDTRIP_FAIL")
