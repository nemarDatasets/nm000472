[![DOI](https://img.shields.io/badge/DOI-10.82901%2Fnemar.nm000472-blue)](https://doi.org/10.82901/nemar.nm000472)

# Example stereo-EEG recording with one clinically confirmed seizure (Zenodo 22150094)

Stereo-EEG (sEEG) recording of one patient during clinical intracranial monitoring, released on Zenodo by
"AI-is-out-there" as "example sEEG (stereo EEG is a specific subtype of intracranial EEG + MRI) with clinical confirmed one seizure" (doi:10.5281/zenodo.22150094, CC-BY-4.0,
version 2 (MRI added)). The record's whole description reads: "'sEEG-HFOs-8.edf': seizure 17:27:11,
clinicaly confirmed 17:27:14 / MRI T1 and T2". No article, protocol, ethics statement or demographic information is
given in the record; this README does not add any.

## Content
- `sub-01/ieeg/sub-01_task-seizure_ieeg.edf`: the released EDF+ file, byte-identical copy (sha256
  16d2be49ccae55379531286de94768f6321ff3ac71860b1c65e41b2a96b502c5). 102 signals (100 sEEG contacts labelled `EEG <electrode><contact>`,
  2 marker inputs `MKR1+`/`MKR2+`) at 256 Hz, 3.92 h (14095 one-second records), physical range
  +-3200 uV, 16-bit. The EDF identification fields are the depositor's placeholders (`patient_2_3_17 X 01-JAN-1999 name`, `Startdate 01-JAN-2019 01.00.00 X X X`),
  so the recording date/time in the header is not the real one.
- Electrodes (12): R, FP, FD, PM, CC, SA, PA, CR', CC', PM', SA', PA'. Names ending in `'` follow the usual sEEG convention for the
  contralateral hemisphere, but the record does not state which hemisphere is which, so no hemisphere is asserted here.
- `events.tsv`: the 3 EDF+ annotations of the file, verbatim (Russian, Windows-1251 in the file, UTF-8 here) with an
  English rendering added by the converter. They mark the clinical onset and a seizure with bilateral tonic-clonic
  evolution ("приступ + БТКП").
  - 10392.734 s: где тут начало?
  - 10396.445 s: приступ + БТКП
  - 10399.469 s: клиника
- Not included: `MRI-with-electrodes.zip` (DICOM T1 and T2 with implanted electrodes, 432 slices). It is not
  defaced, so it is withheld from this redistribution; it is available from the Zenodo record. No electrode coordinates
  are therefore released (`electrodes.tsv` has n/a coordinates).
- `sourcedata/zenodo-22150094/`: the Zenodo record metadata (JSON).

## Provenance
Converted by the ieeg-nemar-20261005 campaign (lane o5c-B) with `code/o5c_B_conv_generic12.py`; the EDF data records are
unchanged. `code/conversion_report.json` lists the header fields, the annotation texts and the checksums.
Note: the EDF patient-code and patient-name header fields were set to X on 2026-10-12 as part of header de-identification; no field values are reproduced here.
