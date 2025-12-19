# DNA Engineer Agent

An agentic framework for designing and manipulating DNA sequences using natural language instructions. Built for use with Claude Code or similar LLM-based coding assistants.

## Quick Start

### 1. Environment

This project lives in your WSL compbio workspace:
```
/home/cnguy/projects/dna_engineer_agent/
```

Synced via GitHub between home laptop (WSL) and work laptop.

**Current access method:** VSCode tunnel from work laptop → home laptop WSL

### 2. Prerequisites

```bash
# Python 3.8+ required
python --version

# Install Biopython (in WSL terminal via VSCode)
pip install biopython
```

### 3. Initial Setup (One Time)

```bash
cd /home/cnguy/projects
# If cloning fresh:
git clone <your-repo-url> dna_engineer_agent
# Or if extracting from zip:
unzip dna_engineer_agent.zip
mv dna_engineer_agent /home/cnguy/projects/

cd /home/cnguy/projects/dna_engineer_agent
python scripts/validate_environment.py
```

You should see all checks pass:
```
✓ PASS: dependencies
✓ PASS: file_structure
✓ PASS: json_valid
✓ PASS: biopython_ops

🟢 AGENT READY: All checks passed.
```

### 4. Initialize Claude Code

**Option A: VSCode Terminal (via tunnel)**
Open a terminal in VSCode (connected to home laptop via tunnel):
```bash
cd /home/cnguy/projects/dna_engineer_agent
claude
```

**Option B: Direct SSH to WSL**
```bash
ssh cnguy@<home-laptop>
cd /home/cnguy/projects/dna_engineer_agent
claude
```

### 5. First Message to Claude Code

```
Read AGENT_INSTRUCTIONS.md completely. Then read cis_elements/manifest.json and parts_library/linkers.json. Finally, load the project context from projects/aav_vhh_v1/context.json. Run scripts/validate_environment.py to confirm you're ready.
```

### 6. Test Design Request

```
Here is test_data/AAV9_minimal_test.gb. 
Please scan it for cis-elements and report what you find.
Then outline how you would insert a placeholder VHH at the VR-IV region using Design 2 linkers.
```

---

## GitHub Sync Workflow

After making changes (new parts, updated cis-elements, etc.):

```bash
cd /home/cnguy/projects/dna_engineer_agent
git add -A
git commit -m "Added new VHH clones to parts library"
git push
```

On work laptop (if working locally there too):
```bash
git pull
```

---

## Directory Structure

```
/home/cnguy/projects/dna_engineer_agent/
├── AGENT_INSTRUCTIONS.md     # Core agent behavior specification
├── README.md                 # This file
├── .gitignore                # Git ignore rules
├── cis_elements/
│   └── manifest.json         # Risk-tagged regulatory elements
├── parts_library/
│   ├── schema.json           # Part definition format
│   └── linkers.json          # Linker designs (D1, D2, D7, etc.)
├── projects/
│   └── aav_vhh_v1/
│       └── context.json      # AAV VHH project-specific rules
├── scripts/
│   └── validate_environment.py  # Environment validation
└── test_data/
    └── AAV9_minimal_test.gb  # Sample GenBank file for testing
```

---

## Core Concepts

### Code-First Philosophy
The agent NEVER manipulates DNA strings directly in conversation. All operations happen through Python/Biopython scripts. This ensures:
- **Reproducibility**: Scripts can be re-run with modified parameters
- **Precision**: No copy-paste errors on large sequences
- **Traceability**: The script IS the documentation

### Mandatory Annotation
Every modification is annotated with:
```
[Action]: <what> | [Reason]: <why> | [Risk]: <level> | [Date]: <when>
```

### Risk Levels
| Level | Description | Agent Behavior |
|-------|-------------|----------------|
| CRITICAL | ITRs, packaging signals | Hard block |
| HIGH | Overlapping ORFs (AAP) | Warn + require confirmation |
| MEDIUM | Promoters, regulatory | Warn but proceed |
| LOW | Common restriction sites | Annotate only |

### Never Overwrite
Output uses versioned filenames: `Parent_v02_Description.gb`

---

## Extending the Agent

### Adding Parts
Create a new JSON file in `parts_library/`:

```json
{
  "library_name": "My VHH Clones",
  "parts": [
    {
      "id": "VHH_CLONE_A",
      "name": "Anti-GFP VHH Clone A",
      "category": "vhh",
      "sequence_type": "protein",
      "sequence_aa": "QVQLVESGGGLVQPGGSLRLSCAASGFTFS...",
      "sequence_dna": null,
      "preferred_codon_table": "Homo sapiens"
    }
  ]
}
```

### Adding Cis-Elements
Edit `cis_elements/manifest.json` to add project-specific regulatory elements:

```json
{
  "id": "MY_CUSTOM_PROMOTER",
  "name": "Custom Promoter",
  "type": "Promoter",
  "detection_method": "sequence_exact",
  "sequence_dna": "TATAATGGGC...",
  "risk_level": "MEDIUM",
  "action_on_overlap": "WARN"
}
```

### Creating New Projects
Create a new directory under `projects/` with a `context.json`:

```
projects/
└── my_new_project/
    └── context.json
```

Then tell the agent: `Load project context from projects/my_new_project/`

---

## Example Workflows

### 1. Scan a Plasmid for Hidden Elements
```
Load test_data/AAV9_minimal_test.gb and scan for all cis-elements.
Report any CRITICAL or HIGH risk features.
```

### 2. Insert a VHH at VR-IV
```
Using the AAV9 file, insert VHH_CLONE_A at the VR-IV region.
Use Design 2 linkers (asymmetric).
Check that AAP reading frame is preserved.
```

### 3. Generate a Library
```
Create three variants of AAV9 with Clone A, B, and C at VR-IV.
Use Design 2 linkers for all.
Output as AAV9_CloneA_D2_v01.gb, AAV9_CloneB_D2_v01.gb, AAV9_CloneC_D2_v01.gb.
```

---

## Troubleshooting

### "Module Bio not found"
```bash
pip install biopython
```

### "JSON validation failed"
Check for trailing commas or missing quotes in your JSON files.

### Agent doesn't find an element
- Verify the sequence is in `cis_elements/manifest.json`
- Check if detection_method matches your use case
- For protein motifs, ensure the reading frame is correct

---

## License

This framework is for research use. No warranty implied.

---

## Version History

- **v1.0.0** (2024-12-16): Initial architecture with AAV VHH v1 project context
