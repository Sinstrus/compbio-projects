const { html2pptx } = require('../../html2pptx/index.cjs');
const PptxGenJS = require('pptxgenjs');
const fs = require('fs');
const path = require('path');

async function buildPresentation() {
  console.log('Building AAV-VHH Display Platform POC presentation...\n');

  // Create new presentation
  const pres = new PptxGenJS();
  pres.layout = 'LAYOUT_16x9';
  pres.title = 'AAV-VHH Display Platform - Phase I POC';
  pres.author = 'Chris Nguyen';

  // Process each slide
  for (let i = 1; i <= 13; i++) {
    const slideNum = String(i).padStart(2, '0');
    const slidePath = path.join(__dirname, `slide${slideNum}.html`);

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

  console.log(`\nTotal slides: 13`);
  console.log('Saving presentation...\n');

  try {
    const outputPath = path.join(__dirname, 'AAV_VHH_Display_POC.pptx');
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
