# AAV Transfer Plasmid Assembly Report v2.0
**Date:** 1766069946.0543737
**Status:** Assembly Complete

## Summary

Successfully assembled two AAV transfer plasmids expressing AAV9 VP1/VP2/VP3 capsid proteins with 6 pre-engineered silent restriction sites and 2 additional junction sites.

## Generated Files

### ssAAV

- **File:** `pGS-ssAAV-EF1A-VP1-rBG_v02.gb`
- **Size:** 8613 bp
- **Cassette:** ~3613 bp
- **Status:** ✅ READY

### scAAV

- **File:** `pGS-scAAV-EF1A-VP1-rBG_v02.gb`
- **Size:** 8576 bp
- **Cassette:** ~3576 bp
- **Status:** ⚠️  WARNING

## Verification Results

| Check | ssAAV | scAAV |
|-------|-------|-------|
| VP1 CDS | ✅ | ✅ |
| VP2 CDS | ✅ | ✅ |
| VP3 CDS | ✅ | ✅ |
| AAP CDS | ✅ | ✅ |
| VR Regions | 9/9 | 9/9 |
| Restriction Sites | 6 | 6 |
| Within Limit | ✅ | ⚠️ |

## Features

### Engineered Restriction Sites

**Internal Sites (from VP1 source):**
1. SmaI - Silent mutation CTT→CTC (L→L)
2. BbvCI - Silent mutation TCC→TCA (S→S)
3. AgeI - Silent mutation ACG→ACC (T→T)
4. BsrGI - Silent mutation GTC→GTA (V→V)
5. BmtI/NheI - Silent mutation GCC→GCT (A→A)
6. BstZ17I - Silent mutation GGA→GGT (G→G)

**Junction Sites (for modular cloning):**
7. AflII (upstream) - Between promoter and VP1
8. NotI (downstream) - Between VP1 stop and 3'UTR

### Expression Cassette Components

```
[EF1α Promoter] → [AflII] → [Kozak+ATG] → [VP1/VP2/VP3/AAP]
    → [NotI] → [mini 3'UTR] → [rBG polyA]
```

## Recommendations

- ✅ **ssAAV construct** is ready for production
- ⚠️  **scAAV construct** may be oversized - test packaging efficiency

## Next Steps

1. Transfect into HEK293T cells
2. Validate VP1/VP2/VP3 expression by Western blot
3. Test capsid assembly
4. Verify restriction site uniqueness by digest

