// =============================================================================
// MathJax 3 Configuration — TrajGRN-Bench
// Works with pymdownx.arithmatex (generic: true)
// =============================================================================
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex",
  },
};
