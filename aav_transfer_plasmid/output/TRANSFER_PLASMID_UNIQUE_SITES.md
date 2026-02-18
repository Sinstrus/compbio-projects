# Transfer Plasmid Unique Restriction Sites - Complete Analysis

**Date:** 2025-12-18
**Approach:** Constrained ONLY by transfer plasmid uniqueness (no RepCap constraints)

---

## Executive Summary

When we remove the RepCap uniqueness constraint and only require uniqueness in the transfer plasmids (ssAAV and scAAV), we have **MANY excellent options** for Regions 1-3:

- **Region 1:** 21 perfectly unique enzymes (all require just 1-2 mutations)
- **Region 2:** 18 perfectly unique enzymes
- **Region 3:** 19 perfectly unique enzymes

All sites require only 1-2 silent mutations and are guaranteed unique in both transfer plasmid contexts.

---

## Region 1: Rep68-stop to VP2-start

**VP1 positions:** 50-410 (0-indexed, 360 bp)

### Perfectly Unique Options (21 total)

#### 1-Mutation Sites (14 options) ⭐ EASIEST

| Enzyme | Type | Applications | Availability |
|--------|------|--------------|--------------|
| **BsaI** | Type IIS | **Golden Gate**, scarless cloning | Excellent |
| **BsmBI** | Type IIS | **Golden Gate**, scarless cloning | Excellent |
| **Esp3I** | Type IIS | **Golden Gate** (same as BsmBI) | Excellent |
| **BspEI** | Type IIS | Directional cloning | Very Good |
| **AvrII** | Standard | General cloning, generates 5' overhang | Excellent |
| **HindIII** | Standard | **Very common**, general cloning | Excellent |
| **KpnI** | Standard | **Very common**, general cloning | Excellent |
| **SalI** | Standard | **Very common**, general cloning | Excellent |
| **XbaI** | Standard | **Very common**, general cloning | Excellent |
| **MluI** | Standard | General cloning, methylation-sensitive | Very Good |
| **NcoI** | Standard | Translation start, maintains ATG | Very Good |
| **HpaI** | Standard | Blunt ends | Good |
| **SnaBI** | Standard | Blunt ends | Good |
| **BstBI** | Standard | Generates 5' overhang | Good |

#### 2-Mutation Sites (7 options)

| Enzyme | Type | Applications | Availability |
|--------|------|--------------|--------------|
| **ClaI** | Standard | General cloning | Very Good |
| **EcoRV** | Standard | Blunt ends, common | Excellent |
| **NsiI** | Standard | General cloning | Good |
| **SphI** | Standard | General cloning | Good |
| **AscI** | 8-bp cutter | Rare cutter, very specific | Good |
| **PmeI** | 8-bp cutter | Rare cutter, blunt ends | Good |
| **SbfI** | 8-bp cutter | Rare cutter | Good |

### Recommended Choices

**Best for Golden Gate:** BsaI, BsmBI, or Esp3I (Type IIS)
**Best for standard cloning:** HindIII, KpnI, SalI, or XbaI (very common)
**Best for rare cutting:** AscI, PmeI, or SbfI (8-bp cutters)

---

## Region 2: VP2-AAP Intergenic

**VP1 positions:** 400-550 (0-indexed, 150 bp)

### Perfectly Unique Options (18 total)

#### 1-Mutation Sites (7 options) ⭐ EASIEST

| Enzyme | Type | Applications | Availability |
|--------|------|--------------|--------------|
| **BsaI** | Type IIS | **Golden Gate**, scarless cloning | Excellent |
| **BsmBI** | Type IIS | **Golden Gate**, scarless cloning | Excellent |
| **Esp3I** | Type IIS | **Golden Gate** (same as BsmBI) | Excellent |
| **BspEI** | Type IIS | Directional cloning | Very Good |
| **KpnI** | Standard | **Very common**, general cloning | Excellent |
| **SalI** | Standard | **Very common**, general cloning | Excellent |
| **XbaI** | Standard | **Very common**, general cloning | Excellent |

#### 2-Mutation Sites (11 options)

| Enzyme | Type | Applications | Availability |
|--------|------|--------------|--------------|
| **AvrII** | Standard | General cloning, 5' overhang | Excellent |
| **ClaI** | Standard | General cloning | Very Good |
| **EcoRV** | Standard | Blunt ends, common | Excellent |
| **HindIII** | Standard | **Very common** | Excellent |
| **HpaI** | Standard | Blunt ends | Good |
| **MluI** | Standard | Methylation-sensitive | Very Good |
| **NcoI** | Standard | Translation start | Very Good |
| **NsiI** | Standard | General cloning | Good |
| **SphI** | Standard | General cloning | Good |
| **BstBI** | Standard | 5' overhang | Good |
| **SbfI** | 8-bp cutter | Rare cutter | Good |

### Recommended Choices

**Best for Golden Gate:** BsaI, BsmBI, or Esp3I
**Best for standard cloning:** KpnI, SalI, or XbaI (1 mutation, very common)
**If 2 mutations OK:** HindIII or EcoRV (extremely common)

---

## Region 3: AAP-stop to VR4

**VP1 positions:** 1120-1350 (0-indexed, 230 bp)

### Perfectly Unique Options (19 total)

#### 1-Mutation Sites (13 options) ⭐ EASIEST

| Enzyme | Type | Applications | Availability |
|--------|------|--------------|--------------|
| **BsaI** | Type IIS | **Golden Gate**, scarless cloning | Excellent |
| **BsmBI** | Type IIS | **Golden Gate**, scarless cloning | Excellent |
| **Esp3I** | Type IIS | **Golden Gate** (same as BsmBI) | Excellent |
| **ClaI** | Standard | General cloning | Very Good |
| **EcoRV** | Standard | Blunt ends, **very common** | Excellent |
| **HindIII** | Standard | **Very common**, general cloning | Excellent |
| **HpaI** | Standard | Blunt ends | Good |
| **KpnI** | Standard | **Very common**, general cloning | Excellent |
| **MluI** | Standard | Methylation-sensitive | Very Good |
| **NcoI** | Standard | Translation start | Very Good |
| **NsiI** | Standard | General cloning | Good |
| **SalI** | Standard | **Very common**, general cloning | Excellent |
| **SnaBI** | Standard | Blunt ends | Good |

#### 2-Mutation Sites (6 options)

| Enzyme | Type | Applications | Availability |
|--------|------|--------------|--------------|
| **AvrII** | Standard | General cloning, 5' overhang | Excellent |
| **BspEI** | Type IIS | Directional cloning | Very Good |
| **BstBI** | Standard | 5' overhang | Good |
| **SphI** | Standard | General cloning | Good |
| **XbaI** | Standard | **Very common** | Excellent |
| **PacI** | 8-bp cutter | Rare cutter | Good |

### Recommended Choices

**Best for Golden Gate:** BsaI, BsmBI, or Esp3I (1 mutation)
**Best for standard cloning:** EcoRV, HindIII, KpnI, or SalI (1 mutation, very common)
**Best for rare cutting:** PacI (8-bp cutter, 2 mutations)

---

## Cross-Region Analysis

### Enzymes Available in ALL Three Regions

These enzymes can be engineered in any of the three regions, giving you maximum flexibility:

| Enzyme | Type | Mutations Needed | Best Choice |
|--------|------|------------------|-------------|
| **BsaI** | Type IIS | 1 in all regions | **Excellent - Golden Gate** |
| **BsmBI** | Type IIS | 1 in all regions | **Excellent - Golden Gate** |
| **Esp3I** | Type IIS | 1 in Regions 2-3, 1 in R1 | **Excellent - Golden Gate** |
| **KpnI** | Standard | 1 in all regions | **Excellent - Very common** |
| **SalI** | Standard | 1 in all regions | **Excellent - Very common** |
| **MluI** | Standard | 1 in all regions | Very Good |
| **NcoI** | Standard | 1 in all regions | Very Good |
| **HpaI** | Standard | 1 in all regions | Good |
| **ClaI** | Standard | 1-2 mutations | Very Good |
| **NsiI** | Standard | 1-2 mutations | Good |

### Enzymes in Two Regions

| Enzyme | Available In | Notes |
|--------|-------------|-------|
| **XbaI** | R1, R2 (1 mut); R3 (2 mut) | Very common, excellent |
| **HindIII** | R1 (1 mut); R2, R3 (1-2 mut) | Very common, excellent |
| **BspEI** | R1, R2 (1 mut); R3 (2 mut) | Type IIS |
| **AvrII** | R1 (1 mut); R2, R3 (2 mut) | Standard |
| **EcoRV** | R1, R2 (2 mut); R3 (1 mut) | Very common |

---

## Recommended Site Combinations

### Option A: Golden Gate Assembly Strategy ⭐

**Focus:** Type IIS enzymes for scarless cloning

| Region | Enzyme | Mutations | Rationale |
|--------|--------|-----------|-----------|
| 1 | **BsaI** | 1 | Type IIS, standard Golden Gate |
| 2 | **BsmBI** | 1 | Type IIS, different from BsaI |
| 3 | **Esp3I** | 1 | Type IIS, compatible |

**Advantages:**
- All Type IIS (cut outside recognition)
- Scarless cloning capability
- Modular assembly workflows
- Only 1 mutation per site

---

### Option B: Standard Cloning - Maximum Availability

**Focus:** Very common enzymes available in every lab

| Region | Enzyme | Mutations | Rationale |
|--------|--------|-----------|-----------|
| 1 | **HindIII** | 1 | Extremely common |
| 2 | **KpnI** | 1 | Extremely common |
| 3 | **SalI** | 1 | Extremely common |

**Advantages:**
- All enzymes ubiquitous in molecular biology labs
- Cheap, reliable, well-characterized
- Only 1 mutation per site
- Generate sticky ends for traditional cloning

---

### Option C: Mixed Strategy - Best of Both

**Focus:** Combine Golden Gate with standard enzymes

| Region | Enzyme | Mutations | Rationale |
|--------|--------|-----------|-----------|
| 1 | **BsaI** | 1 | Golden Gate in early region |
| 2 | **XbaI** | 1 | Standard, very common |
| 3 | **EcoRV** | 1 | Blunt ends, very common |

**Advantages:**
- Golden Gate option in Region 1
- Standard enzymes for Regions 2-3
- Maximum flexibility
- Only 1 mutation per site

---

### Option D: All-BsaI Strategy (if compatible)

**Simplest:** Use BsaI in all three regions

| Region | Enzyme | Mutations | Rationale |
|--------|--------|-----------|-----------|
| 1 | **BsaI** | 1 | Type IIS |
| 2 | **BsaI** | 1 | Type IIS |
| 3 | **BsaI** | 1 | Type IIS |

**Note:** Need to verify that the three BsaI sites are in different sequence contexts to avoid unintended cutting. BsaI is Type IIS and cuts outside its recognition site, so positional context matters for Golden Gate.

---

## Implementation Priority

### Highest Priority (1 mutation, very common):

1. **HindIII, KpnI, SalI, XbaI** - Universal availability
2. **BsaI, BsmBI, Esp3I** - Golden Gate capability
3. **EcoRV** - Blunt ends, very common

### Medium Priority (1 mutation, specialized):

4. **MluI, NcoI, ClaI** - Good enzymes, slightly less common
5. **BspEI** - Type IIS, directional cloning

### Lower Priority (2 mutations):

6. **8-bp cutters (AscI, PmeI, PacI, SbfI)** - Rare cutters for special applications

---

## Next Steps

1. **Choose strategy:** Golden Gate (Type IIS), Standard (common enzymes), or Mixed
2. **Select specific enzymes** for each of the 3 regions
3. **Design silent mutations** for each chosen site
4. **Verify mutations are silent** in VP1 reading frame
5. **Update RepCap plasmid** with new restriction sites
6. **Generate v03 transfer plasmids**
7. **Validate uniqueness** in final constructs

---

## Key Advantages of This Approach

✅ **Maximum flexibility** - 50+ total options across 3 regions
✅ **All require ≤2 mutations** - Easy to engineer
✅ **All perfectly unique** - Guaranteed no conflicts in transfer plasmids
✅ **Type IIS options** - Golden Gate assembly capability
✅ **Common enzymes** - Universal lab availability
✅ **8-bp cutters** - Rare cutting for special needs

---

**Analysis Complete** ✓
