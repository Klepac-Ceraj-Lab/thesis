# Getting the data

From the root directory of this project, run


```sh
# metadata
curl -o diabimmune/data/DIABIMMUNE_Karelia_metadata.RData https://diabimmune.broadinstitute.org/diabimmune/uploads/attachments/87/DIABIMMUNE_Karelia_metadata.RData

# amplicon data
curl -o diabimmune/data/DIABIMMUNE_Karelia_16S_data.RData https://diabimmune.broadinstitute.org/diabimmune/uploads/attachments/68/diabimmune_karelia_16s_data.rdata

# mgx data
curl -o diabimmune/data/DIABIMMUNE_karelia_metaphlan_table.txt https://diabimmune.broadinstitute.org/diabimmune/uploads/attachments/70/diabimmune_karelia_metaphlan_table.txt
```