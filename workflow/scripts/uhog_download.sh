# download uhgg cripsr spacers host taxonomy
wget -O /home/carsonjm/resources/bacteria_db/uhog_genome-all_metadata.tsv http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-oral/v1.0/genomes-all_metadata.tsv

# iterate through metadata and extract genome file paths
awk '{FS="\t"} (NR>1) {print $20}' /home/carsonjm/resources/bacteria_db/uhog_genome-all_metadata.tsv > /home/carsonjm/resources/bacteria_db/uhog_genome-ftp_list
sed -i 's/ftp:/http:/' /home/carsonjm/resources/bacteria_db/uhog_genome-ftp_list

# remove trailing character
sed -i 's/\r$//' /home/carsonjm/resources/bacteria_db/uhog_genome-ftp_list

# download uhgg genomes in parallel
cat /home/carsonjm/resources/bacteria_db/uhog_genome-ftp_list | xargs -n 1 -P 5 wget -P /home/carsonjm/resources/bacteria_db/bacteria -q -w 5 -T 5

# create file indicating download is complete
touch /home/carsonjm/resources/bacteria_db/bacteria/bacteria_db_download_complete