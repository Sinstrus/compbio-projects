const { html2pptx } = require('../../html2pptx/index.cjs');
const PptxGenJS = require('pptxgenjs');
const fs = require('fs');
const path = require('path');

async function buildPresentation() {
  console.log('Building AAV-VHH Display Platform presentation (Voyager theme)...\n');

  const pres = new PptxGenJS();
  pres.layout = 'LAYOUT_16x9';
  pres.title = 'AAV-VHH Display Platform - Voyager Theme';
  pres.author = 'Chris Nguyen';
  pres.subject = 'AAV Nanobody Display for Targeted Gene Therapy';

  const slides = [
    "slide01_title.html",
    "slide02_problem.html",
    "slide03_solution.html",
    "slide04_section_design.html",
    "slide05_insertion_sites.html",
    "slide06_plasmid_table.html",
    "slide07_section_production.html",
    "slide08_strategies.html",
    "slide09_exp_groups.html",
    "slide10_section_validation.html",
    "slide11_pipeline.html",
    "slide12_timeline.html",
    "slide13_summary.html",
  ];

  for (let i = 0; i < slides.length; i++) {
    const slideNum = String(i + 1).padStart(2, '0');
    const slidePath = path.join(__dirname, slides[i]);

    console.log(`  Processing slide ${slideNum}...`);

    try {
      await html2pptx(slidePath, pres, {
        tmpDir: fs.realpathSync(path.join(__dirname, 'assets')),
      });
      console.log(`    ✓ Slide ${slideNum} added`);
    } catch (error) {
      console.error(`    ✗ Error on slide ${slideNum}:`, error.message);
      throw error;
    }
  }

  console.log(`\nTotal slides: ${slides.length}`);
  console.log('Saving presentation...\n');

  try {
    const outputPath = path.join(__dirname, 'AAV_VHH_Display_Voyager.pptx');
    await pres.writeFile({ fileName: outputPath });

    const stats = fs.statSync(outputPath);
    console.log(`SUCCESS: Presentation saved to ${outputPath}`);
    console.log(`File size: ${(stats.size / 1024).toFixed(1)} KB`);
  } catch (error) {
    console.error('ERROR saving file:', error.message);
    process.exit(1);
  }
}

buildPresentation();
