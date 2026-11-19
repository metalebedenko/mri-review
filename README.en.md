# An assistant for reviewing MRI scans

[По-русски](README.md)

`mri-review` is a skill for an AI agent that explains a radiologist’s report in plain language. The written report alone is enough for that. When reviewing images, the agent can point to the series and frame behind an observation. For scans taken on different dates, it compares matching series and notes questions to discuss with the radiologist.

This is a research prototype. Its clinical accuracy has not been validated. A radiologist should review any medical conclusions drawn from the images.

## What an explanation looks like

*A fictional line from a report*

> “Assessment of this series is partially limited by motion artefact. No prior studies were available.”

*An example answer to “What does that mean?”*

> Movement blurred some of the images. The radiologist did not have an earlier scan for comparison, so this sentence cannot tell you what has changed over time. One question to ask is whether the area of interest is clear enough on these images.

## Using the skill

Place the `mri-review` folder in the user skills directory of a compatible AI agent. Start with [SKILL.md](SKILL.md). If your agent supports calling skills by name, use `$mri-review`.

Tell the agent what you have and what you want to understand. For example, “I only have the written report. Please explain it in plain language and list questions I can ask the radiologist.”

## Where measurements come from

The skill includes a separate calculator for distances. A person selects landmarks on a single frame and supplies the image geometry. The calculator converts coordinates into millimetres. It supports lengths, distances to an infinite line and projections along a patient axis. A DICOM viewer and an anatomical recognition model are outside the package.

The [synthetic example](references/synthetic-landmarks.json) requires Python 3.10 or later. Run this command from the project directory. It uses the standard library and works offline.

```sh
python3 scripts/measure_landmarks.py references/synthetic-landmarks.json
```

The expected distance to the line is 3 mm, the length is 5 mm, and the signed projection along the chosen axis is −4 mm (magnitude 4 mm). Points from different frames require a different calculation. [Data format and geometry](references/measurements.md) · [Source code](scripts/measure_landmarks.py) · [Arithmetic checks](tests/)

## Data and privacy

The repository contains synthetic data only. Keep actual scans, annotations and results elsewhere. The calculator runs locally and retains the original `input_record` fields without anonymising them. Before sharing material with an AI agent, check where that agent sends data. [Data handling](SECURITY.md)

By default, results go to the screen. The `--output` option creates a new file on a filesystem that supports hard links. It refuses to overwrite an existing file.

## Further reading

- [Reviewing and comparing studies](references/workflow.md) · [Anatomical profiles](references/anatomy-profiles.md)
- [Specialised MRI series](references/specialized-mri.md) · [Writing up observations and questions](references/reporting-and-followup.md)
- [MIT licence for the original code and text](LICENSE) · [Third-party material and rights](THIRD_PARTY_NOTICES.md)
