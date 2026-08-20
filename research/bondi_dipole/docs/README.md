# Working notes — Bondi dipole

These are the **working documents** behind the campaign: how the matter was
made to behave, what was tried and rejected, and the running journal.  None of
them is cited by the paper, and none is part of the published results pack —
that is why they live here and not in `results/bondi-dipole-runaway/`, which is
kept to the artefacts the article actually needs.

Read them if you are rebuilding this work or trying to understand *why* the
setup looks the way it does.  Read the paper if you want the results.

| file | what it is |
|---|---|
| [`FINDINGS.md`](FINDINGS.md) | the physics as it stood during the campaign: results, tests against the point-mass reference, caveats, open questions. Largely superseded by the article, kept for the reasoning trail |
| [`DEBUGGING.md`](DEBUGGING.md) | how the matter was made to stop dispersing — two root causes and the falsifications that found them. Not in the paper, and the most useful document here if you are reproducing the setup |
| [`Debug.md`](Debug.md) | the long-form scrutiny log: every objection raised against the result and what was done about it |
| [`GPU.md`](GPU.md) | what still has to be run on a GPU and why: the queue for the next campaign, what each item unblocks, what blocks it, and what costs how much. Also what *not* to spend device time on |
| [`campaign_journal.md`](campaign_journal.md) | the campaign journal, verbatim and unedited |

## Where everything else lives

* **The paper** — [`../bondi_dipole.tex`](../bondi_dipole.tex)
* **The figures and the script that builds them** — [`../figures/`](../figures/),
  [`../make_article_figures.py`](../make_article_figures.py)
* **The published results pack** — [`../../../results/bondi-dipole-runaway/`](../../../results/bondi-dipole-runaway/),
  whose own [`README.md`](../../../results/bondi-dipole-runaway/README.md) and
  [`campaign/README.md`](../../../results/bondi-dipole-runaway/campaign/README.md)
  document every run and where it is used in the article.
