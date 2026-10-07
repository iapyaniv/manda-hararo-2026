"""Estimate the AGU counted words and publication units for a LaTeX source.

AGU counts the abstract, body text, section headings, in-text citations,
figure captions, acknowledgments, and appendices; each equation counts as
one word. It excludes front matter, the Plain Language Summary, table
contents, Open Research, references, and supporting information.
Publication units = words/500 + figures + tables; GRL allows 12.
https://www.agu.org/publications/authors/journals/publication-fees

Requires TeXcount (texcount.pl): https://app.uio.no/ifi/texcount/

Usage:
    python scripts/agu_word_count.py --texcount /path/to/texcount.pl main.tex

This is a source estimate, not a submission count: citation commands are
not expanded (the bibliography is not read), and unwritten additions are
not included.
"""
import argparse
import re
import subprocess
import tempfile
from pathlib import Path


FIELDS = (
    'Words in text', 'Words in headers',
    'Words outside text (captions, etc.)',
    'Number of math inlines', 'Number of math displayed',
)

# Sections AGU does not count, each ignored up to the next section command.
EXCLUDED_SECTIONS = (r'Plain Language Summary', r'Open Research Section')


def counting_view(source):
    # Original %TC directives (e.g. abstract-only markers) can hide the
    # main text entirely, so replace them with explicit boundaries.
    text = re.sub(r'^[ \t]*%TC:.*\n', '', source, flags=re.M)
    start = re.search(r'^\\begin\{abstract\}', text, re.M)
    end = re.search(r'^\\end\{document\}', text, re.M)
    if not start or not end or start.start() >= end.start():
        raise SystemExit('Expected an abstract followed by \\end{document}')
    text = text[start.start():end.start()]
    for heading in EXCLUDED_SECTIONS:
        pattern = (r'(\\section\*?\{' + heading + r'\}.*?)(?=\\section[*{])')
        text, matches = re.subn(
            pattern, r'%TC:ignore\n\1%TC:endignore\n', text, flags=re.S)
        if matches > 1:
            raise SystemExit(f'More than one "{heading}" section')
    return ('\\documentclass{article}\n\\begin{document}\n'
            + text + '\\end{document}\n')


def run_texcount(texcount, tex_path):
    result = subprocess.run(
        ['perl', str(texcount), '-inc', str(tex_path)],
        check=True, capture_output=True, text=True)
    if result.stderr:
        raise SystemExit(result.stderr)
    labels = FIELDS + ('Number of floats/tables/figures',)
    values = {
        key: int(re.search(
            r'^' + re.escape(key) + r': (\d+)$', result.stdout, re.M
        ).group(1))
        for key in labels
    }
    sections = []
    for match in re.finditer(
            r'^\s+(\d+)\+(\d+)\+(\d+) \((\d+)/(\d+)/(\d+)/(\d+)\) (.+)$',
            result.stdout, re.M):
        total = sum(int(match.group(i)) for i in (1, 2, 3, 6, 7))
        label = match.group(8).split('}\\label{')[0].replace('_top_', 'Abstract')
        sections.append((label, total))
    return values, sections


def main():
    parser = argparse.ArgumentParser(
        description='Estimate AGU counted words for a LaTeX manuscript.')
    parser.add_argument('--texcount', type=Path, required=True,
                        help='path to texcount.pl')
    parser.add_argument('source', type=Path, help='manuscript .tex file')
    args = parser.parse_args()
    view = counting_view(args.source.read_text())
    with tempfile.TemporaryDirectory(prefix='agu-word-count-') as work:
        path = Path(work) / args.source.name
        path.write_text(view)
        values, sections = run_texcount(args.texcount, path)
    words = sum(values[key] for key in FIELDS)
    figures = values['Number of floats/tables/figures']
    print(f'Source: {args.source}')
    print('Scope: abstract through end of document; Plain Language Summary '
          'and Open Research excluded; citations not expanded.')
    for label, total in sections:
        print(f'  {total:>6}  {label}')
    print(f'Counted words (text + headings + captions + equations): {words}')
    print(f'Floats (figures and tables, including placeholders): {figures}')
    print(f'Publication units (words/500 + floats): '
          f'{words / 500 + figures:.1f} of 12 allowed (GRL)')


if __name__ == '__main__':
    main()
