# Contributing

Thanks for helping improve this course. Students and interns are encouraged to contribute — spotting a confusing explanation is just as valuable as adding new material.

## Quick fixes (typos, unclear wording)

Open the file on GitHub, click the pencil icon, make the change and propose it as a pull request.

## Questions or larger suggestions

Open an [issue](https://github.com/kanji709/genetic-evaluation-101/issues) describing:

- which lesson (and script or section of `tutorial.qmd`) it relates to
- what was unclear or wrong
- (optional) how you would fix it

## Changing the course

1. Fork the repo and create a branch, e.g. `git checkout -b lesson-03-reml-note`
2. Follow "Editing the course" in the README:
   - change functions in `R/functions.R`, then run `Rscript tools/sync_functions.R`
   - edit text and exercises directly in `tutorial.qmd`
   - if you change the simulation, delete `tutorial_cache/` and run `Rscript 00_simulate_data.R`
3. Check that `quarto render tutorial.qmd` succeeds and the lesson scripts run in order from the repository root
4. Open a pull request with a short description of what you changed

## Data

Use **simulated or openly licensed data only**. Never commit real animal, farm or client records.
