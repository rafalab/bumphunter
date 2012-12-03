(TeX-add-style-hook "bumphunter"
 (lambda ()
    (LaTeX-add-bibliographies)
    (TeX-add-symbols
     '("bold" 1)
     '("software" 1)
     '("Rpackage" 1)
     '("Rcode" 1)
     "R")
    (TeX-run-style-hooks
     "natbib"
     "numbers"
     "fullpage"
     "hyperref"
     "color"
     "times"
     "latex2e"
     "art12"
     "article"
     "12pt")))

