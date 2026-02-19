const { html2pptx } = require('../../html2pptx/index.cjs');
const PptxGenJS = require('pptxgenjs');
const fs = require('fs');
const path = require('path');

async function buildPresentation() {
  console.log('Building N-Terminal Intron Retention presentation (Voyager theme)...\n');

  // Create new presentation
  const pres = new PptxGenJS();
  pres.layout = 'LAYOUT_16x9';
  pres.title = 'N-Terminal Intron Retention: VP1 and VP2 VHH Display';
  pres.author = 'Payload Group';
  pres.subject = 'Extending Intron Retention to N-Terminal Positions';
  pres.company = 'Voyager Therapeutics';

  const slides = [
    "slide01-title.html",
    "slide02-executive-summary.html",
    "slide03-rationale.html",
    "slide04-cassette-architecture.html",
    "slide05-avd132-vp2.html",
    "slide06-avd133-vp1.html",
    "slide07-splicing-topology.html",
    "slide08-risk-assessment.html",
    "slide09-verification.html",
    "slide10-conclusions.html"
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
    const outputPath = path.join(__dirname, 'NTerm_Intron_Retention_VP1_VP2_Voyager.pptx');
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
