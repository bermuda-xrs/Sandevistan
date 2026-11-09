import katex from "katex";

/** Render a TeX string with the configured KaTeX instance. Display only; no arithmetic. */
export function renderFormula(tex: string, displayMode = false): string {
  return katex.renderToString(tex, {
    throwOnError: false,
    displayMode,
    output: "html",
  });
}
