# Collaboration Plan

## Team Members
- George: Led data loading and initial BEM solver implementation.
- Milton: Developed airfoil polar interpolation routines and performance curve plotting.
- Vaso: Managed documentation and testing framework.

## Roles and Responsibilities
Each member owned specific modules and supported one another:
- *George*: implemented data_loader.py and the blade geometry parser.
- *Milton*: created airfoil.py and integrated it with the solver.
- *Vaso*: organized the repository layout and wrote examples/main.py.

## Communication
- Primary channel: WhatsApp for daily status updates and quick questions.

## Workflow
1. Tasks were documented as GitHub Issues and assigned to team members.
2. Branches followed the convention feature/<task-name> and <task-name>/develop.
3. Each feature branch was merged into main only after a peer review and approval.
4. Before merging, unit tests were added or updated, and code style checks (pylint, black) were run locally.

## Meetings
- Weekly meetings on Mondays at 17:00 for progress reviews and task assignments.
- Ad‑hoc in‑person sessions at the library after lectures.

## Code Review and Quality Assurance
- Pull requests required at least one approving review.
- CI pipeline on GitHub Actions ran pytest --cov=src tests/ and pylint src/ on every PR.
- Achieved over 80% test coverage and maintained a pylint score above 8.0.

## Tooling
- *Version Control*: GitHub for code hosting and issue tracking.
- *Communication*: WhatsApp.
- *Continuous Integration*: GitHub Actions for testing and style checks.

## Acknowledgments
All team members contributed collaboratively, provided constructive feedback, and ensured the codebase met the course requirements and quality standards.