# Create or amend the documentation

You should draft documentation for the wiki for any major new features that need explaining, or update the existing docs for any changes.

The docs are currently on the branch `enhancement/documentation` but will soon move to the `main/develop` branch. For now you can check out this branch and edit and push to it without review (assuming you have access).

Once you have the relevant branch you can view the docs locally in your browser by installing mkdocs:

```pip install mkdocs```

and then running

```mkdocs serve```

in the folder in which the index `mkdocs.yml` lives (the base folder).

It should say something like `Serving on http://127.0.0.1:8000/`, and you can then open this address in your browser to view the docs as they will appear online once you push them after editing.