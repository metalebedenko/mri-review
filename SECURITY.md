# Security and data handling

This repository is a research example. It contains no patient studies, medical images, reports, identifiers, or access credentials. Examples and automated tests must remain synthetic.

## Data boundary

- The Python calculator reads an explicitly supplied JSON file and writes to stdout or a new output file. It performs no network requests.
- Output includes the complete `input_record`, including optional and unknown metadata. It is **not a de-identification tool**. Input and output may contain sensitive information if you supply it.
- Keep real inputs, outputs, DICOM files, screenshots, reports, credentials, and environment files outside this repository. Ignore rules are a guardrail, not an anonymization or security audit.
- An AI agent following the skill may use external models or tools. Their data handling is separate from this local calculator. Do not assume local processing or send sensitive material without an appropriate authorization and data-handling arrangement.
- The arithmetic does not establish anatomical correctness, diagnostic accuracy, or clinical suitability. Do not use its result alone for medical decisions.

## Reporting

For a reproducible problem, create an issue using a minimal **synthetic** example. Do not include a patient case, screenshot, identifier, local private path, token, or real medical file. Do not put exploit details involving other systems or credentials in a public issue. A maintainer can reproduce numerical or parsing problems without private materials.

There is no guaranteed response time or supported clinical deployment.
