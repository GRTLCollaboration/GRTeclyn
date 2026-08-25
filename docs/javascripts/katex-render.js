document.addEventListener('DOMContentLoaded', function () {
  function renderMath() {
    renderMathInElement(document.body, {
      delimiters: [
        { left: '\\(', right: '\\)', display: false },
        { left: '\\[', right: '\\]', display: true },
      ],
      throwOnError: false,
    });
  }

  renderMath();

  // Re-run on mkdocs theme navigation (readthedocs uses AJAX page loads)
  var observer = new MutationObserver(function () {
    renderMath();
  });
  observer.observe(document.body, { childList: true, subtree: true });
});
