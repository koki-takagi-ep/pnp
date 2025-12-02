# Repository Guidelines

## Project Structure & Module Organization
- Root docs: writing and reference files. Key files: `survey.md` (PNP数値解法サーベイ), `implicit_newton_design.md` (設計メモ), `AGENTS.md` (本ガイド)。
- `survey/`: 研究PDFと対応する`.bib`（例: `ding_2019_jcp.pdf` / `ding_2019_jcp.bib`）と参照集 `references.bib` を格納。
- Keep new PDFs and their BibTeX pairs inside `survey/` using `familyname_year_journalabbr` naming (e.g., `takagi_2024_jap.pdf` / `.bib`).

## Build, Test, and Development Commands
- No build system here; content is Markdown + Bib + PDFs. Validate structure with quick checks:
  - List survey files: `ls survey`
  - Search text: `rg "pattern" survey`

## Coding Style & Naming Conventions
- Markdown: use `#`, `##` for headings; keep sections concise.
- BibTeX keys and filenames: `familyname_year_journalabbr` (short journal abbr), e.g., `liu_2023_jsc`. No `abstract` fields; include `doi` when available.
- PDFs live beside their `.bib` in `survey/`. Avoid spaces in filenames.

## Testing Guidelines
- No automated tests. For content edits, manually open Markdown for formatting and check BibTeX syntax (commas, braces, required fields).

## Commit & Pull Request Guidelines
- Keep commits focused and descriptive (e.g., `add ding_2019_jcp bib/pdf`, `update survey summary`).
- In PRs: summarize changes (files touched, new references), note naming compliance, and mention if any files were moved.

## Agent Notes
- Do not delete or overwrite existing PDFs/Bib without explicit reason. If structure changes (e.g., new subfolders), update `survey.md` and this guide accordingly.
- When adding references, prefer open-access sources; if download fails, document the attempt and leave placeholders rather than guessing content.
