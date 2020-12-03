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

## Fastq files

Fastq files that have been through kneaddata are also available,
but are too big to be put into a git repo.

There are [instructions here](https://diabimmune.broadinstitute.org/diabimmune/three-country-cohort/resources/metagenomic-sequence-data)
for downloading everything, but they didn't work for me.
I ran

```sh
wget -r -np -nd https://pubs.broadinstitute.org/diabimmune/data/16
```

but this just created a file called `16`,
which appeared to be an html file with a bunch of urls.
This is a bit hacky, but it worked.

```
mv 16 urls.txt
sed -i 's/href/\n/g' 16 # replaces `href` with newlines
sed -Ei "s/^='(.+)' download.+/\1/" urls.txt # replaces everything other than the url
cat urls.txt | xargs -n1 -P12 wget
```

Current path for these files is on the lab server `hopper`,
at `/augusta/students/danielle/diabimmune`.