const { html2pptx } = require('../../html2pptx/index.cjs');
const PptxGenJS = require('pptxgenjs');
const fs = require('fs');
const path = require('path');

async function buildPresentation() {
  console.log('Building TRACER-Nano VR4 Splicing Validation presentation (Voyager theme)...\n');

  // Create new presentation
  const pres = new PptxGenJS();
  pres.layout = 'LAYOUT_16x9';
  pres.title = 'TRACER-Nano VR4 Splicing Design Validation';
  pres.author = 'Payload Group';
  pres.subject = 'VR4 Loop-Based VHH Display';
  pres.company = 'Voyager Therapeutics';

  const slides = [
    "slide01-title.html",
    "slide02-executive-summary.html",
    "slide03-vr4-concept.html",
    "slide04-config-comparison.html",
    "slide05-splice-donor.html",
    "slide06-splice-acceptor.html",
    "slide07-frame-preservation.html",
    "slide08-stoichiometry.html",
    "slide09-cryptic-sites.html",
    "slide10-risk-assessment.html",
    "slide11-validation-plan.html",
    "slide12-conclusions.html"
  ];

  // Process each slide
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
    const outputPath = path.join(__dirname, 'TRACER_Nano_VR4_Splicing_Slides.pptx');
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
