## Summary

I've successfully created a skill proposal for your project that addresses genome quality assessment for downstream processing. Here's what I accomplished:

### Created Skill Proposal
- **Skill Name**: `genome-quality-assessor`
- **Proposal ID**: `genome-quality-assessor-20260818-0ab71f66e6` (currently pending)
- **Purpose**: Assess genome quality for downstream processing by providing three key metrics: genome completeness, contamination level, and confidence scores

### Key Features Addressing Your Requirements
1. **Completeness Score**: Using BUSCO to measure genome completeness as percentage of expected single-copy orthologs found
2. **Contamination Assessment**: Identification and quantification of likely contaminant sequences
3. **Confidence Metric**: Reliability score for downstream predictions based on genome quality

### Integration with Your Use Cases
The skill is designed to work with your four benchmark applications:
1. **Rhizobia/Paraburkholderia**: Compare ANI vs mobile nod/nif functions with quality awareness
2. **Burkholderia/Paraburkholderia**: Distinguish environmental vs clinical lineages reliably
3. **Arbuscular Mycorrhizal Fungi**: Assess fragmentation effects on taxonomy and gene presence
4. **Metagenome-Assembled Genomes**: Develop completeness-aware interpretation for host-microbe interactions

### Connection to Existing Skills
Designed to integrate with:
- `busco-assessor`: For detailed completeness assessment
- `genome-compare`: For ANI-based relatedness analysis with quality awareness
- `metagenomics`: For functional capacity assessment

### Next Steps
The skill proposal is currently pending review. Once approved, it can be implemented and will be available via:
```
clawbio run genome-quality-assessor --demo
```

I've also created a detailed markdown proposal file (`GENOME_QUALITY_ASSESSOR_PROPOSAL.md`) for your review before any implementation begins.

This skill directly answers your project goal of determining "how good is a genome for downstream processing" by providing the three essential numbers you requested, with explicit quality-aware interpretation to prevent false conclusions from fragmented or contaminated assemblies.