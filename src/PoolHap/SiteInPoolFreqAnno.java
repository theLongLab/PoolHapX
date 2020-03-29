package PoolHap;

import java.io.BufferedReader;
import java.io.FileReader;

public class SiteInPoolFreqAnno {

    public int num_sites;
    public int num_pools;
    public String[] pool_IDs;
    public LocusAnnotation[] loci_annotations;
    public double[][] inpool_freqs; // [this.num_sites][this.num_pools]

    /**
     *  Gold standard variant positions file format below:
     *  (First line is the headerheader; the rest are variant annotations (column #1) and
     *  frequencies.)
     *
     *  e.g.:
     *  Var_ID\t0\t1\t2\t3\t4
     *  0;329;329;0:1\t0.8054936896807721\t0.8954512105649303\t1.0\t0.7920611798980335\t0.807418
     *  0;796;796;0:1\t0.1968\t0.10485592315901815\t0.0\t0.20496397117694154\t0.1938584779706275
     *  0;833;833;0:1\t0.20240320427236316\t0.2144571885836223\t0.0952253934382502\t0.09605122732123
     *
     *  Gold standard variant frequency annotation object constructor.
     *
     *  @param input_file (required) gold standard variant positions file path string.
     */
    public SiteInPoolFreqAnno(String input_file) {
        try {
            /*
             *  Read the file to determine the lengths of related fields.
             *  TODO: [ReconEP]:: this section should be a helper I think.
             */
            BufferedReader br = new BufferedReader(new FileReader(input_file));
            String line = br.readLine();

            // Skip metadata lines (starting with "#"), if any.
            // TODO: [ReconEP]:: extract to static utils method.
            while (line.startsWith("#")) {
                line = br.readLine();
            }

            // Determine number of pools from length of header line (minus annotation column).
            String[] header = line.split("\t");
            this.num_pools = header.length - 1;
            line = br.readLine();

            // Determine number of sites from number of data lines.
            // TODO: [ReconEP]:: mentioned in GC file; occurs multiple times, extract to global
            // variable?
            while (line != null) {
                this.num_sites++;
                line = br.readLine();
            }

            br.close();


            /*
             *  Initialize array variables based on previously determined lenghts.
             */
            this.inpool_freqs = new double[this.num_sites][this.num_pools];
            this.pool_IDs = new String[this.num_pools];
            for (int p = 0; p < this.num_pools; p++) {
                this.pool_IDs[p] = header[p + 1];
            }

            this.loci_annotations = new LocusAnnotation[this.num_sites];


            /*
             *  Read the file a second time to fill arrays.
             *  TODO: [ReconEP]:: I think this part can also be a separate helper.
             */
            br = new BufferedReader(new FileReader(input_file));
            line = br.readLine();

            // Skip metadata lines.
            // TODO: [ReconEP]:: as per above, extract to static utils method.
            while (line.startsWith("#")) {
                line = br.readLine();
            }

            line = br.readLine(); // skip the header
            int locus_index = 0;

            // Parse data lines.
            while (line != null) {
                String[] tmp = line.split("\t");
                this.loci_annotations[locus_index] = new LocusAnnotation(tmp[0]);

                // TODO: [LEFTOVER]
                // System.out.println("ding\t" + this.loci_annotations[0].alleles_coding.size());

                for (int p = 0; p < this.num_pools; p++) {
                    this.inpool_freqs[locus_index][p] = Double.parseDouble(tmp[p + 1]);
                }

                locus_index++;
                line = br.readLine();
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

    }
    
    
    public SiteInPoolFreqAnno(String input_file, int start, int end, int pool_index ) {
        try {
            /*
             *  Read the file to determine the lengths of related fields.
             *  TODO: [ReconEP]:: this section should be a helper I think.
             */
            BufferedReader br = new BufferedReader(new FileReader(input_file));
            String line = br.readLine();

            // Skip metadata lines (starting with "#"), if any.
            // TODO: [ReconEP]:: extract to static utils method.
            while (line.startsWith("#")) {
                line = br.readLine();
            }

            // Determine number of pools from length of header line (minus annotation column).
            String[] header = line.split("\t");
            this.num_pools =  header.length - 1;
            line = br.readLine();

            // Determine number of sites from number of data lines.
            // TODO: [ReconEP]:: mentioned in GC file; occurs multiple times, extract to global
            // variable?
            int pos_index =0;
            while (line != null) {
            	if ((pos_index>= start) && (pos_index<= end )) {
                	this.num_sites++;
                	line = br.readLine();
            	}
            	pos_index++;
            }

            br.close();


            /*
             *  Initialize array variables based on previously determined lenghts.
             */
            this.inpool_freqs = new double[this.num_sites][this.num_pools];
            this.pool_IDs = new String[this.num_pools];
            for (int p = 0; p < this.num_pools; p++) {
            	this.pool_IDs[p] = header[p + 1];
            }

            this.loci_annotations = new LocusAnnotation[this.num_sites];


            /*
             *  Read the file a second time to fill arrays.
             *  TODO: [ReconEP]:: I think this part can also be a separate helper.
             */
            br = new BufferedReader(new FileReader(input_file));
            line = br.readLine();

            // Skip metadata lines.
            // TODO: [ReconEP]:: as per above, extract to static utils method.
            while (line.startsWith("#")) {
                line = br.readLine();
            }

            line = br.readLine(); // skip the header
            int locus_index = 0;

            // Parse data lines.
            while (line != null) {
                String[] tmp = line.split("\t");
                if ((locus_index>= start) && (locus_index<= end )) {
	                this.loci_annotations[locus_index-start] = new LocusAnnotation(tmp[0]);
	
	                // TODO: [LEFTOVER]
	                // System.out.println("ding\t" + this.loci_annotations[0].alleles_coding.size());
	
	                for (int p = 0; p < this.num_pools; p++) {
	                    this.inpool_freqs[locus_index-start][p] = 
	                    		Double.parseDouble(tmp[p + 1]);
	                }
                }
                locus_index++;
                line = br.readLine();
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

    }


    /**
     *  Constructor: copy by reference -- no clone!
     *  Gold standard variant frequency annotation copy constructor.
     *  Copy by reference -- no clone!
     *
     *  @param num_sites
     *  @param num_pools
     *  @param pool_IDs
     *  @param loci_annotations
     *  @param inpool_freqs
     */
    public SiteInPoolFreqAnno(
        int num_sites,
        int num_pools,
        String[] pool_IDs,
        LocusAnnotation[] loci_annotations,
        double[][] inpool_freqs) {

        this.num_sites = num_sites;
        this.num_pools = num_pools;
        this.pool_IDs = pool_IDs;
        this.loci_annotations = loci_annotations;
        this.inpool_freqs = inpool_freqs;
    }


    /**
     *  Extract the user-specified subset of loci by index and return a new object. To save memory,
     *  here the annotations and frequencies are copied by reference.
     *
     *  @param start_index
     *  @param end_index
     */
    public SiteInPoolFreqAnno subset(int start_index, int end_index) {
        int sub_num_sites = end_index - start_index + 1;
        LocusAnnotation[] sub_loci_annotations = new LocusAnnotation[sub_num_sites];
        double[][] sub_inpool_freq = new double[sub_num_sites][];
        for (int i = start_index; i <= end_index; i++) {
            sub_loci_annotations[i - start_index] = this.loci_annotations[i];
            sub_inpool_freq[i - start_index] = this.inpool_freqs[i];
        }

        return new SiteInPoolFreqAnno(
            sub_num_sites,
            this.num_pools,
            this.pool_IDs,
            sub_loci_annotations, sub_inpool_freq);

    }
}
