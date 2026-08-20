# Pull requests

Here we describe the rough procedure that one should follow for making changes to the code. It is not set in stone, but strongly advised as best practice.

No one is too new or junior to add code to the repository, you don't need permission, but if you are unsure contact Katy or Miren for informal guidance. It may take a little time for it to be reviewed, but please always consider adding your code if you have made a nice new feature, and consider it *an absolute moral obligation* to report bugs (by creating a [new issue](https://github.com/GRTeclyn/GRTeclyn/issues/new/choose)). Although contributing back your changes might not naively seem like the absolute best use of your time, remember that if everyone contributes back their features and modifications, this saves time for others when they need them for a project later on; this is the power of open source coding!

## Step 1: Testing

Thoroughly test your code and ensure as far as practically possible that it does not introduce any regressions. If you are pushing to a branch in this repo (see [below](#step-2-a-fresh-copy) on how to do this), then the tests are automatically run with multiple compilers every time you push new changes, although this is not sufficient on its own to ensure correctness of your new feature and that it does not introduce regressions.  For large new features or substantial modifications, ideally one article should have been published using it before it is ready to be merged into this repository. Equally, if you have new code from a published article, please do make it public! If it was useful for you it will certainly be useful for others. GRTeclyn is open source for a reason.

## Step 2: A fresh copy

If your fork of the GRTeclyn repository is messy (as many of ours are) and has significant divergences from the `main` branch in this repository, it might be easiest to start with a new branch based off this `main` branch and copy in your changes incrementally to this new branch. Make sure you have fetched/pulled the latest commits from the `main` branch on GitHub before starting. When pushing to GitHub, please push to this repository rather than your own fork.

> **Note**
> Contact [Katy](https://github.com/KAClough/) or
> [Miren](https://github.com/mirenradia) if you don't have the necessary
> permissions to push to this repository.

If you have a local clone of your fork, you can add the main public "upstream" repository as a new remote repository using the command
```bash
git remote add upstream https://github.com/GRTLCollaboration/GRTeclyn.git
```
You can replace `upstream` with whatever name you prefer e.g. `public` but we shall stick with `upstream` below. Then, make sure you fetch the latest commits from this new remote using the command
```bash
git fetch upstream
```
Finally, to create a new branch based off upstream's `main`, you can do the two commands
```bash
git checkout upstream/main
git branch enhancement/my_new_enhancement
```
Finally, once you have added and committed your changes to this new branch, you can push it to `upstream` by using the command
```bash
git push -u upstream enhancement/my_new_enhancement
```
This local branch will be set to track a branch of the same name in the upstream repository (rather than on your fork).

Alternatively, if you don't want to have to deal with multiple remote repositories, you can clone this repository to a separate directory to your fork.

It is probably useful to make your changes in several staged commits, so that you can look back through the process of what you did. For example, if you plan to add a new file and then amend it, first add the file, then commit, then amend it, then commit again. This way the commit history will show the changes you made. If you only commit after editing, the whole thing will appear as a new file, which isn't very helpful.

When writing commit messages, try to follow [this style guide](https://cbea.ms/git-commit/). In particular, rather than
writing one long message using a command such as
```bash
git commit -m "My really really really long single-line commit message that goes on for far too long"
```
do
```bash
git commit
```
which should open you into a text editor where you can separate out the short title of the commit (ideally <= 50 chars) from a more detailed description in the body (this is not always necessary for small changes).

> **Note**
> The default editor is usually `vim` but if
> you prefer a different editor you can change it using the command
> ```bash
> git config --global core.editor "<my preferred text editor command>"
> ```

Small commits for minor changes (e.g. fixing a single typo) should be avoided. These can be rebased or squashed into one commit using an
interactive rebase:
```bash
git rebase -i
```
but if you have the code already developed, you should be able to avoid them by adding things carefully in a logical order.

Having said the advice above, don't be too fussy; it's better to get the (tidy, readable) code in with a messy commit history than not add it at all.

## Step 3: Compile and debug

Unless your code was already up-to-date with the main branch, you may now have compilation issues that you need to fix. You should add a new test for any major changes to confirm that the feature still behaves as expected and later changes don't introduce regressions. Looking at the existing tests will give you ideas on how to set these up, but ask for help if you aren't sure.

## Step 4: Tidy and comment the code

All new code will be held to high standards of readability. You should follow the design principles in this wiki, and the [**coding conventions**](coding_conventions.md). In particular you will need to ensure code is formatted using Clang Format.

## Step 5: Create or amend the documentation

You should also draft documentation for the wiki for any major new features that need explaining, or update the existing docs for any changes.

See [**Updating the documentation**](updating_docs.md) for more info on how to do this.

## Step 6: Create a pull request (PR)

Open a [new Pull request (PR)](https://github.com/GRTLCollaboration/GRTeclyn/pulls) in this repository, requesting to merge code from your new branch into the `main` branch.

This is a place for your new code to be reviewed by other members of the collaboration and improved prior to inclusion. The format is mostly self-explanatory, but will highlight the changes to the code you have made, and show whether the automatically run unit tests have passed. More details on GitHub pull requests are given here:
[Pull Requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests)

The main thing you need to do is write some text which explains what your amendment does, why it is needed, and how you have confirmed that it works as required. You should also explain the reasoning behind any unusual approaches you may have taken.

You can make the PR
[draft](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests#draft-pull-requests)
if you don't want anyone to review it yet but might appreciate constructive comments.

## Step 7 : Review

The code needs to be reviewed by someone else in the collaboration, and strictly it should not be someone from the same institution as you or who has collaborated on the same project that produced the code. It needs to be someone critical who can look at it with fresh eyes and really challenge whether it is good quality code. The person you want to review it should be tagged in the PR, if in doubt tag both @mirenradia and @KAClough and we will do it or delegate.


> TOP TIP: Don't take the review personally! Use it as a learning experience and really engage in the discussions to make the code as good as it can be. (Respectful) creative conflict is your friend.