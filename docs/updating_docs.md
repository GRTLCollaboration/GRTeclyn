# Create or amend the documentation

You should draft documentation for the wiki for any major new features that need explaining, or update the existing docs for any changes.

The docs are currently on the branch `enhancement/documentation` but will soon move to the `main/develop` branch. For now you can check out this branch and edit and push to it without review (assuming you have access).

Once you have the relevant branch you can view the docs locally in your browser. First, create and activate a Python virtual environment using your preferred tool, for example:

```
python3 -m venv .venv
source .venv/bin/activate
```

Then install the required Python packages from `docs/requirements.txt`:

```
pip install -r docs/requirements.txt
```

and then running

```mkdocs serve```

in the folder in which the index `mkdocs.yml` lives (the base folder).

It should say something like `Serving on http://127.0.0.1:8000/`, and you can then open this address in your browser to view the docs as they will appear online once you push them after editing.

## Writing mathematical expressions

The docs support rendering of LaTeX math expressions using [KaTeX](https://katex.org/) via the [`pymdownx.arithmatex`](https://facelessuser.github.io/pymdown-extensions/extensions/arithmatex/) extension. Use standard LaTeX delimiters:

- **Inline math**: `$...$`, e.g. `$\chi$` renders as $\chi$
- **Display (block) math**: `$$...$$` on its own line, e.g.:

  ```
  $$
  H^2 = \frac{8 \pi G \rho}{3}
  $$
  ```

Do **not** use the GitHub-style backtick syntax (`$`...`$`) as it will not render correctly with mkdocs.