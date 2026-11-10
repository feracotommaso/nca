NCA REVISIONS: COLOR SUPPORT FOR QUARTO PDF AND DOCX

Files:
- nca_paper_R1_colors_enabled.qmd  : your updated manuscript, with the color settings enabled
- nca_revision_colors.lua         : Pandoc filter for blue/red inline revision spans
- nca_revision_reference.docx     : Word reference template with Revision Blue and Revision Red character styles

Place ALL THREE files in the same folder. Render the colors_enabled QMD
with your existing Quarto project and apaquarto extension:

  quarto render nca_paper_R1_colors_enabled.qmd --to apaquarto-pdf
  quarto render nca_paper_R1_colors_enabled.qmd --to apaquarto-docx

The manuscript's revision markup is currently BLUE:
  [changed text]{.revision style="color: blue"}

Change blue to red in a span if you prefer red for that passage.

IMPORTANT:
- Keep your existing project materials (e.g. nca_refs.bib, data files,
  scripts, and apaquarto extension). They are not included in this ZIP.
- The Word reference file is based on Pandoc's default reference.docx.
  If your APA extension requires a particular Word reference template,
  it may need its own styles copied into the provided file, or vice versa.
- The only manuscript changes in the provided QMD are its active filter
  and the DOCX reference-doc setting. Your article text is unchanged.
- A PDF renderer requires xcolor, which your QMD already imports.
