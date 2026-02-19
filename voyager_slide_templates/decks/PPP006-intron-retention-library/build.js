const { html2pptx } = require('../../html2pptx/index.cjs');
const PptxGenJS = require('pptxgenjs');
const fs = require('fs');
const path = require('path');

async function buildPresentation() {
  console.log('Building VHH3 Intron Retention Library presentation (Voyager theme)...\n');

  // Create new presentation
  const pres = new PptxGenJS();
  pres.layout = 'LAYOUT_16x9';
  pres.title = 'VHH3 Intron Retention Library: Splice Engineering at VR4';
  pres.author = 'Payload Group';
  pres.subject = 'Intron Retention Cassettes for VHH Display on AAV9';
  pres.company = 'Voyager Therapeutics';

  const slides = [
    "slide01-title.html",
    "slide02-executive-summary.html",
    "slide03-concept.html",
    "slide04-cassette-architecture.html",
    "slide05-u1-framework.html",
    "slide06-tiered-donors.html",
    "slide07-hek293t-strategy.html",
    "slide08-avd-library.html",
    "slide09-cryptic-problem.html",
    "slide10-dnachisel-pipeline.html",
    "slide11-before-after.html",
    "slide12-verification.html",
    "slide13-next-steps.html",
    "slide14-conclusions.html"
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
    const outputPath = path.join(__dirname, 'VHH3_Intron_Retention_Library.pptx');
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
