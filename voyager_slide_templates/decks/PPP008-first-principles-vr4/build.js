const { html2pptx } = require('../../html2pptx/index.cjs');
const PptxGenJS = require('pptxgenjs');
const fs = require('fs');
const path = require('path');

async function buildPresentation() {
  console.log('Building First-Principles VR4 Splice Donor Library presentation...\n');

  const pres = new PptxGenJS();
  pres.layout = 'LAYOUT_16x9';
  pres.title = 'First-Principles Splice Donor Library for VR4 Intron Retention';
  pres.author = 'Payload Group';
  pres.subject = 'AVD100-131: 32 Synthetic Splice Donor Constructs';
  pres.company = 'Voyager Therapeutics';

  const slides = [
    "slide01-title.html",
    "slide02-executive-summary.html",
    "slide03-motivation.html",
    "slide04-donor-panel.html",
    "slide05-cassette-architecture.html",
    "slide06-design-comparison.html",
    "slide07-protein-products.html",
    "slide08-build-results.html",
    "slide09-qc-verification.html",
    "slide10-next-steps.html"
  ];

  const tmpDir = fs.realpathSync(path.join(__dirname, 'assets'));

  for (let i = 0; i < slides.length; i++) {
    const slideNum = String(i + 1).padStart(2, '0');
    const slidePath = path.join(__dirname, slides[i]);

    console.log(`  Processing slide ${slideNum}...`);

    try {
      await html2pptx(slidePath, pres, { tmpDir });
      console.log(`    ✓ Slide ${slideNum} added`);
    } catch (error) {
      console.error(`    ✗ Error on slide ${slideNum}:`, error.message);
      throw error;
    }
  }

  console.log(`\nTotal slides: ${slides.length}`);
  console.log('Saving presentation...\n');

  try {
    const outputPath = path.join(__dirname, 'First_Principles_VR4_Splice_Donor_Library.pptx');
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
