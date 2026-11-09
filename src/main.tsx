import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import "katex/dist/katex.min.css";
import { renderFormula } from "./formula/renderFormula.ts";
import App from "./App.tsx";
import "./index.css";

// KaTeX is loaded with the page. M0 does not draw a formula yet.
void renderFormula;

createRoot(document.getElementById("root")!).render(
  <StrictMode>
    <App />
  </StrictMode>,
);
