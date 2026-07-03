// =============================================================================
// Extra JavaScript for TrajGRN-Bench documentation
// =============================================================================

// Make Plotly charts responsive
window.addEventListener("resize", function () {
  if (typeof Plotly !== "undefined") {
    var charts = document.querySelectorAll(".js-plotly-plot, .plotly-graph-div");
    charts.forEach(function (chart) {
      Plotly.Plots.resize(chart);
    });
  }
});
