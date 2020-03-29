package PoolHap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import PoolHap.LocusAnnotation;

public class HapConfig {
    /**
     *  @author  Quan Long. Oct 12, 2018
     *
     *  Recording the configuration of the haplotypes, both in-pool and global.
     *
     *  The objects can be a small region or a whole chromosome.
     *  The object of this class will be input or output parameters of other relaying algorithms.
     *
     *  Corresponding to the data structure in this class, there are two files:
     *  (1) global_hap_file:
     * 	    (1.1) Each row represents a locus; each column represents a haplotype.
     * 	    (1.2) The first row is the header of the IDs of all haps (usually are indexes);
     *            the second row is the global frequency (or NaN if unknown).
     * 	    (1.3) The first column is the IDs of all loci.
     * 	          (1.3.1) chr_index;start_loc;end_loc;alleles (the alleles are separated by ":").
     *                    chr_index starts with zero.
     *            (1.3.2) NOTE: if it is an indel, start_loc=end_loc; Here start_loc!=end_loc only
     *                    if it is a region, instead of a primitive locus.
     * 	    (1.4) Use '\t' to separate columns.
     *  (2) local_hap_file: local_haplotypes
     * 	    (2.1) Each row represents the in-pool frequencies of all haplotypes (many of them can be
     *            zero)
     * 	          the first row is the header of all hap IDs (in the global file)
     * 	          the first element is the pool ID
     * 	          the rest elements are the frequencies of this haplotype in all pools: zero
     *            indicates the absence.
     *      (2.2) Each column represents a haplotype's frequencies in all pools.
     *            the first column is the IDs of all pools.
     * 	    (2.3) Use '\t' to separate columns.
     *      (2.4) In general, it is recommended that the order of haplotypes in this file is the
     *            same as the order of the haplotypes in the global hap file. But if it is not the
     *            case, the program can handle it.
     *      (Note that this matrix representation is convenient for analysis of relationship between
     *      pools, e.g., transmission or evolution.)
     *
     *  Note that the information of how to encode alleles to numbers will be generated on-the-fly
     *  in the future, we may add a function to read them from a file.
     */

    /*
     *  Redundant variables for convenience.
     */
    public int num_global_hap; // number of global haplotypes in this region
    public int num_loci; // number of loci in this region
    public int num_pools; // number of pools under study

    // Map the pool IDs to their indices in this.hap_IDs.
    public HashMap<String, Integer> hapID2index;

    /*
     *  Main structures.
     */
    public String[] hap_IDs; // # of global_hap

    // # of global_hap x # of loci: haplotypes that ever show up in the total population. String
    // coded.
    public String[][] global_haps_string;

    // # of global_hap x # of loci: haplotypes that ever show up in the total population. Floating
    // number coded.
    public double[][] global_haps;
    public double[] global_haps_freq; // # of global_hap
    public String[] pool_IDs; // # of pools
    public double[][] in_pool_haps_freq; // # of global_hap x # of pools
    public LocusAnnotation[] locusInfo; // # of loci; note that a locus can be a SNP or a region
    public double[][] inpool_site_freqs; // # of loci x # of pools; Added by Quan Dec. 2018.
    boolean Using1asDiagonal;

    /*
     *  Summary statistics.
     */
    // The 'average' allele at each position i.e. average haplotype Average frequency of alternate
    // allele.
    public double[] mu;
    public double[][] sigma; // the variance-covariance matrix for each position

    // The log-likellihood of this HapConfig explaining the observed variant data.
    public double logL;

    // Estimated number of individuals per pool (for now, assumed same per pool but easily
    // extended). Used to calculate logL.
    public int est_ind_pool;

    /*
     *  Link to the solving method
     */
    RegionEMSolver regional_EM_solver;

    // Constructors

    /**
     *  Constructor, forming the object by the variables in the memory.
     *  This constructor is useful in the divide & conquer algorithm.
     *
     *  @param global_haps_string global haplotype matrix
     *  @param global_haps_freq frequencies of global haplotypes double array
     *  @param in_pool_haps_freq frequencies of haplotypes in-pool matrix
     *  @param inpool_site_freqs frequencies of sites in-pool matrix
     *  @param locusInput annotated loci array
     *  @param num_pools number of pools
     *  @param hap_IDs haplotype ID array
     *  @param pool_IDs pool ID array
     *  @param est_ind_pool // TODO: [Question]:: what is this?
     */
    public HapConfig(
        String[][] global_haps_string,
        double[] global_haps_freq,
        double[][] in_pool_haps_freq,
        double[][] inpool_site_freqs,
        LocusAnnotation[] locusInput,
        int num_pools,
        String[] hap_IDs,
        String[] pool_IDs,
        int est_ind_pool) {

        // Set known variables
        this.num_global_hap = global_haps_freq.length;
        this.num_loci = locusInput.length;
        this.global_haps_freq = global_haps_freq.clone();
        this.global_haps_string = global_haps_string.clone();
        for (int k = 0; k < this.num_global_hap; k++) {
            this.global_haps_string[k]=global_haps_string[k].clone();
            // TODO: [LEFTOVER] Quan to Michael: I have moved the setID to a function, if this is what you meant.
            // System.out.println(this.hap_IDs[k]);
        }
        if (hap_IDs != null) { // *** Changed from this.hap_IDs (global) to hap_IDs (parameter)
            this.hap_IDs = hap_IDs.clone();
        } else { // if no IDs assigned, use its haplotype 0/1 string for the moment. 
        	this.set_binary_hap_IDs_using_hapString();        	
        }
        if (in_pool_haps_freq != null) { // *** Matched format as above
            this.in_pool_haps_freq = in_pool_haps_freq.clone();
            for (int j = 0; j < this.num_global_hap; j++) {
                for (int k = 0; k < this.num_pools; k++) {
                    this.in_pool_haps_freq[j][k] = in_pool_haps_freq[j][k];
                }
            }
        } else { // if no intra-pool frequencies provided, start at 0
            this.in_pool_haps_freq = new double[this.num_global_hap][this.num_pools];
        }
        if (inpool_site_freqs != null) {
            this.inpool_site_freqs = inpool_site_freqs.clone();
        }

        // Originally, clone(), which is not a deep-clone. But we assume that the locusInfo won't be
        // changed in the algorithm.
        this.locusInfo = locusInput;

        // // Need to specifically copy the mappings over to the new locus storage object.
        // // Shallow copy alone doesn't work.

        this.construct_hapID2index_map();
        this.encoding_haps(); // initialize this.global_haps
        this.num_pools = num_pools;
        if (pool_IDs != null) { // *** matched format as above
            this.pool_IDs = pool_IDs.clone();
        } 
        else {
        	System.out.println("ERROR: pool_IDs == null");
        	System.exit(0);
        }

        this.est_ind_pool = est_ind_pool;
        this.update_sigma_mu_logL();
    }
    
    
    public HapConfig(
            String[][] global_haps_string,
            double[] global_haps_freq,
            double[][] in_pool_haps_freq,
            double[][] inpool_site_freqs,
            LocusAnnotation[] locusInput,
            int num_pools,
            String[] hap_IDs,
            String[] pool_IDs,
            int est_ind_pool,
            boolean ok ) {

            // Set known variables
            this.num_global_hap = global_haps_freq.length;
            this.num_loci = locusInput.length;
            this.global_haps_freq = global_haps_freq.clone();
            this.global_haps_string = global_haps_string.clone();
            for (int k = 0; k < this.num_global_hap; k++) {
                this.global_haps_string[k]=global_haps_string[k].clone();
                // TODO: [LEFTOVER] Quan to Michael: I have moved the setID to a function, if this is what you meant.
                // System.out.println(this.hap_IDs[k]);
            }
            if (hap_IDs != null) { // *** Changed from this.hap_IDs (global) to hap_IDs (parameter)
                this.hap_IDs = hap_IDs.clone();
            } else { // if no IDs assigned, use its haplotype 0/1 string for the moment. 
            	this.set_binary_hap_IDs_using_hapString();        	
            }
            if (in_pool_haps_freq != null) { // *** Matched format as above
                this.in_pool_haps_freq = in_pool_haps_freq.clone();
                for (int j = 0; j < this.num_global_hap; j++) {
                    for (int k = 0; k < this.num_pools; k++) {
                        this.in_pool_haps_freq[j][k] = in_pool_haps_freq[j][k];
                    }
                }
            } else { // if no intra-pool frequencies provided, start at 0
                this.in_pool_haps_freq = new double[this.num_global_hap][this.num_pools];
            }
            if (inpool_site_freqs != null) {
                this.inpool_site_freqs = inpool_site_freqs.clone();
            }
            this.Using1asDiagonal= ok;

            // Originally, clone(), which is not a deep-clone. But we assume that the locusInfo won't be
            // changed in the algorithm.
            this.locusInfo = locusInput;

            // // Need to specifically copy the mappings over to the new locus storage object.
            // // Shallow copy alone doesn't work.

            this.construct_hapID2index_map();
            this.encoding_haps(); // initialize this.global_haps
            this.num_pools = num_pools;
            if (pool_IDs != null) { // *** matched format as above
                this.pool_IDs = pool_IDs.clone();
            } 
            else {
            	System.out.println("ERROR: pool_IDs == null");
            	System.exit(0);
            }

            this.est_ind_pool = est_ind_pool;
            this.update_sigma_mu_logL();
        }
    
    

    /**
     * Use haplotype strings to set a binary ID.
     */
    public void set_binary_hap_IDs_using_hapString() {
    	this.hap_IDs=new String[this.num_global_hap];
    	for(int h=0;h<this.num_global_hap;h++) {
    		String id = "";
            for (int locus = 0; locus < this.num_loci; locus++) {
                id += this.global_haps_string[h][locus];
            }
            this.hap_IDs[h] = id;
    	}
    }
    /**
     *	Combine inpool_haps into a multi-pool HapConfig
     *  Initially by Lauren, modified and commented by Quan 2019-07
     *  Before entering the method, one needs to ensure that only one pool in each inpool_haps[p]
     *  
     *  @param inpool_haps  // renamed from @param final_local_haps in Lauren's version
     */
    public HapConfig(HapConfig[] inpool_haps) {
        this.num_loci = inpool_haps[0].num_loci;
        ArrayList<String> global_ids = new ArrayList<String>();
        ArrayList<ArrayList<String>> local_ids = new ArrayList<ArrayList<String>>();
        for (int p = 0; p < inpool_haps.length; p++) {
            local_ids.add(new ArrayList<String>());
            for (int h = 0; h < inpool_haps[p].num_global_hap; h++) {
                String curr_vc = "";
                for (int l = 0; l < this.num_loci; l++) {
                    curr_vc += inpool_haps[p].global_haps_string[h][l];
                }
                if (!global_ids.contains(curr_vc)) {
                    global_ids.add(curr_vc);
                }
                local_ids.get(p).add(curr_vc);
            }
        }
        this.num_global_hap = global_ids.size();
        //this.hap_IDs = new String[this.num_global_hap];
        this.num_pools = inpool_haps.length;
        this.global_haps_string = new String[this.num_global_hap][this.num_loci];
        this.in_pool_haps_freq = new double[this.num_global_hap][inpool_haps.length];
        double[] tmp_global_freq = new double[this.num_global_hap];
        for (int h = 0; h < this.num_global_hap; h++) {
            String[] hap_var_comp = global_ids.get(h).split("");
            for (int locus = 0; locus < this.num_loci; locus++) {
                this.global_haps_string[h][locus] = hap_var_comp[locus];
            }
            for (int p = 0; p < this.num_pools; p++) {
                if (local_ids.get(p).contains(global_ids.get(h))) {
                    int inpool_index = local_ids.get(p).indexOf(global_ids.get(h));
                    this.in_pool_haps_freq[h][p] =
                    		inpool_haps[p].global_haps_freq[inpool_index];
                    // inpool_haps[p] has only one pool, therefore its global_haps_freq is the inpool frequency
                    tmp_global_freq[h] += inpool_haps[p].global_haps_freq[inpool_index];
                }
            }
            //this.hap_IDs[h] = "h"+h;
        }
        // normalize the global_haps_freq
        this.global_haps_freq = new double[this.num_global_hap];
        for (int h = 0; h < this.num_global_hap; h++) {
            this.global_haps_freq[h] = tmp_global_freq[h] / this.num_pools;
        }
        // set up the rest 
        this.set_binary_hap_IDs_using_hapString();
        this.recode_HapIDs_to_base16();
        this.inpool_site_freqs = inpool_haps[0].inpool_site_freqs.clone();
        this.locusInfo = inpool_haps[0].locusInfo;
        this.construct_hapID2index_map();
        this.encoding_haps();
        this.pool_IDs = new String[this.num_pools];
        for (int p = 0; p < this.num_pools; p++) {
            this.pool_IDs[p] = inpool_haps[p].pool_IDs[0];  // there is only one ID in each pool.
        }
        this.est_ind_pool = 0;
        this.update_sigma_mu_logL();
    }
    
    public HapConfig(String regression_in_file, String hap_out_file, 
    		HashMap<Integer, String> index_var_prefix_dict ) throws IOException {
    	ArrayList<String >  final_haps= new ArrayList<String>();
    	ArrayList<Double >  final_freq= new ArrayList<Double>();
    	File file = new File(regression_in_file); 
    	if (file.exists() ) {
	    	BufferedReader bufferedreader = new BufferedReader(new FileReader(regression_in_file));
	    	String line = "";
	    	double freq_count=0.00001;
	    	boolean is_write= true;
	        while ((line = bufferedreader.readLine()) != null) {
	        	if (!line.substring(0, 1).equals("#")) {
		        	line= line.replace("\n", "").replace("\r", "");
		        	String[] tmp = line.split("\t");
		        	if ((tmp[1].equals("NA") ) || (tmp[1].equals("") )) {
		        		is_write =false;
		        	}
		        	if (is_write ) {
			        	if    (Double.parseDouble(tmp[1]) >0 ) {
			        		final_haps.add(tmp[2]); 
			        		final_freq.add(Double.parseDouble(tmp[1]));
			        		freq_count+=Double.parseDouble(tmp[1]);
			        	}
		        	}
	        	}
	        }
	        bufferedreader.close();
	        
	        if (is_write ) {
//		        freq_count= 1.0;
		        
		        for (int i=0;i< final_freq.size();i++) {
		        	final_freq.set(i, final_freq.get(i)/freq_count);
		        }
		        BufferedWriter bw = new BufferedWriter(new FileWriter(hap_out_file));
		        bw.write("Hap_ID");
		        for (int h = 0; h < final_haps.size(); h++) {
		            bw.write("\th" + Integer.toString(h));
		        }
		
		        bw.write("\nFreq");
		        for (int h = 0; h < final_haps.size(); h++) {
		            bw.write("\t" + final_freq.get(h) );
		        }
		        bw.write("\n");
		        for (int l = 0; l < final_haps.get(0).length(); l++) {
		            bw.write(index_var_prefix_dict.get(l) );
		            for (int h = 0; h < final_haps.size(); h++) {
		                bw.write("\t" + final_haps.get(h).substring(l,l+1));
		            }
		            bw.write("\n");
		        }
		        
		        bw.close();
	        }
    	} 
    	
    }

    /**
     *  Constructor that reads information from files.
     *  This constructor is for the connection between modules (such as AEM, rjMCMC, etc.)
     *
     *  @param global_hap_input_file
     *  @param in_pool_hap_input_file
     */
    public HapConfig(String global_hap_input_file, String in_pool_hap_input_file) {
        try {
            // Parse global_hap_input_file.
            BufferedReader br = new BufferedReader(new FileReader(global_hap_input_file));
            String[] header_ids = br.readLine().split("\t");
            String[] header_freqs = br.readLine().split("\t");
            this.num_global_hap = header_ids.length - 1;
            this.hap_IDs = new String[this.num_global_hap];
            this.global_haps_freq = new double[this.num_global_hap];
            for (int h = 0; h < num_global_hap; h++) {
                this.hap_IDs[h] = header_ids[h + 1];
                this.global_haps_freq[h] = Double.parseDouble(header_freqs[h + 1]);
            }

            String line = br.readLine();
            while (line != null) {
                this.num_loci++;
                line = br.readLine();
            }

            br.close();
            this.global_haps_string = new String[this.num_global_hap][this.num_loci];
            this.locusInfo = new LocusAnnotation[this.num_loci];
            br = new BufferedReader(new FileReader(global_hap_input_file));
            line = br.readLine();
            line = br.readLine(); // skip two headers
            line = br.readLine();
            int loci_index = 0;
            while (line != null) {
                String[] tmp = line.split("\t");

                // The first column is the locus-info.
                this.locusInfo[loci_index] = new LocusAnnotation(tmp[0]);
                for (int h_index = 0; h_index < this.num_global_hap; h_index++) {
                    this.global_haps_string[h_index][loci_index]=tmp[h_index+1];
                }

                // TODO: [LEFTOVER]
                // System.out.println(this.locusInfo[loci_index].alleles_coding.get("1"));

                loci_index++;
                line = br.readLine();
            }

            br.close();
            this.construct_hapID2index_map();
            this.encoding_haps();  // initialize this.global_haps

            // Then read the local-hap file.
            if (in_pool_hap_input_file != null) {            
                br = new BufferedReader(new FileReader(in_pool_hap_input_file));
                String[] header_hap_IDs = br.readLine().split("\t"); // the first row is the hap IDs

                // A *_haps.intra_freq.txt (in-pool frequencies of haplotypes) has been provided.
                if (header_hap_IDs[0].equals("Hap_ID")) {
                    if (this.num_global_hap != header_hap_IDs.length - 1) {
                        System.out.println(
                            "WRONG: this.num_global_hap is not the same between two files!");
                    }
                    line = br.readLine();
                    while (line != null){ // count how many pools
                        this.num_pools++;
                        line = br.readLine();
                    }

                    br.close();
                    this.pool_IDs = new String[this.num_pools];
                    this.in_pool_haps_freq = new double[this.num_global_hap][this.num_pools];
                    br = new BufferedReader(new FileReader(in_pool_hap_input_file));
                    line = br.readLine(); // read the file again, skip the header
                    line = br.readLine();
                    int pool_index=0;
                    while (line != null) {
                        String[] tmp = line.split("\t");
                        this.pool_IDs[pool_index] = tmp[0]; // assign pool IDs
                        for (int h = 1; h < header_hap_IDs.length; h++) {
                            // Based on the hap-ID, find out the hap-index.
                            int h_index = this.hapID2index.get(header_hap_IDs[h]);
                            this.in_pool_haps_freq[h_index][pool_index] =
                                Double.parseDouble(tmp[h]);

                        }
                        pool_index++;
                        line = br.readLine();
                    }
                    br.close();

                // A *_vars.intra_freq.txt (in-pool frequencies of alternate alleles) has been
                // provided.
                } else {
                    this.num_pools = header_hap_IDs.length - 1;
                    this.pool_IDs = new String[this.num_pools];
                    for (int p = 0; p < this.num_pools; p++) {
                        this.pool_IDs[p] = header_hap_IDs[p + 1];
                    }
                    this.inpool_site_freqs = new double[this.num_loci][this.num_pools];
                    int locus_index = 0;
                    line = br.readLine();
                    while(line != null){
                        String tmp[] = line.split("\t");
                        for (int p = 0; p < this.num_pools; p++) {
                            this.inpool_site_freqs[locus_index][p] = Double.parseDouble(tmp[p + 1]);
                        }
                        locus_index++;
                        line = br.readLine();
                    }
                    br.close();
                }
            }

        } catch(Exception e) {
            e.printStackTrace();
        }
    }


    /** This constructor is for haplotypes in GC output files.
     *
     *  @param gc_header
     *  @param var_file
     *  @param num_pools
     *  @throws IOException
     */
//    public HapConfig(
//        String gc_header,
//        String var_file,
//        int num_pools) throws IOException {
//
//        HashSet<String> hap_tracker = new HashSet<String>();
//        ArrayList<HashMap<String, Integer>> freq_tracker =
//            new ArrayList<HashMap<String, Integer>>();
//
//        String[] tmp = null;
//        this.num_pools = num_pools;
//        for (int p = 0; p < this.num_pools; p++) {
//            freq_tracker.add(new HashMap<String, Integer>());
//            //BufferedReader br = new BufferedReader(new FileReader(gc_header +  p + ".in"));
//    			TODO: [Quan] hard-coded file path above!!
//            String line = br.readLine();
//            while (line != null) {
//                tmp = line.split("\\t|-|\\?");
//                String hap = "";
//                for (int l = 0; l < tmp.length - 1; l++) {
//                    hap += tmp[l];
//                }
//                int hap_ct = Integer.parseInt(tmp[tmp.length - 1]);
//                if (!hap_tracker.contains(hap)) {
//                    hap_tracker.add(hap);
//                }
//                if (!freq_tracker.get(p).containsKey(hap)) {
//                    freq_tracker.get(p).put(hap, hap_ct);
//                } else {
//                    int prev_ct = freq_tracker.get(p).get(hap);
//                    freq_tracker.get(p).put(hap, prev_ct + hap_ct);
//                }
//                line = br.readLine();
//            }
//            br.close();
//        }
//
//        this.num_loci = tmp.length - 1;
//        this.num_global_hap = hap_tracker.size();
//        this.global_haps_string = new String[this.num_global_hap][this.num_loci];
//        //this.hap_IDs = new String[this.num_global_hap];
//        int[][] hap_ct = new int[this.num_global_hap][this.num_pools];
//        int hap_index = 0;
//        for (String hap : hap_tracker) {
//            String[] tmp2 = hap.split("");
//            for (int l = 0; l < this.num_loci; l++) {
//                this.global_haps_string[hap_index][l] = tmp2[l];
//            }
//
//            for (int p = 0; p < this.num_pools; p++) {
//                if (freq_tracker.get(p).containsKey(hap)) {
//                    hap_ct[hap_index][p] = freq_tracker.get(p).get(hap);
//
//                } else {
//                    hap_ct[hap_index][p] = 0;
//                }
//
//            }
////            this.hap_IDs[hap_index] = "h"+Integer.toString(hap_index);
//            hap_index++;
//        }
//        this.set_binary_hap_IDs_using_hapString();
//        this.global_haps_freq = new double[this.num_global_hap];
//        this.in_pool_haps_freq = new double[this.num_global_hap][this.num_pools];
//        for (int p = 0; p < this.num_pools; p++) {
//            int pool_size = 0;
//            for (int h = 0; h < this.num_global_hap; h++) {
//                pool_size += hap_ct[h][p];
//            }
//
//            for (int h = 0; h < this.num_global_hap; h++) {
//                this.in_pool_haps_freq[h][p] = (double) hap_ct[h][p] / pool_size;
//            }
//
//        }
//
//        for (int h = 0; h < this.num_global_hap; h++) {
//            double curr_hap_global = 0;
//            for (int p = 0; p < num_pools; p++) {
//                curr_hap_global += this.in_pool_haps_freq[h][p];
//            }
//
//            this.global_haps_freq[h] = curr_hap_global / num_pools;
//        }
//
//        this.locusInfo = new LocusAnnotation[this.num_loci];
//        this.pool_IDs = new String[num_pools];
//        for (int k = 0; k < this.num_pools; k++) {
//            this.pool_IDs[k] = "p"+k;//TODO Quan
//        }
//
//        this.inpool_site_freqs = new double[this.num_loci][this.num_pools];
//        BufferedReader br = new BufferedReader(new FileReader(var_file));
//        String line = br.readLine();
//        line = br.readLine();
//        int locus_index = 0;
//        while (line != null) {
//            tmp = line.split("\t");
//
//            // The first column is the locus-info.
//            this.locusInfo[locus_index] = new LocusAnnotation(tmp[0]);
//            for (int p = 0; p < this.num_pools; p++) {
//                this.inpool_site_freqs[locus_index][p] = Double.parseDouble(tmp[p + 1]);
//            }
//            locus_index++;
//            line = br.readLine();
//        }
//
//        br.close();
//        this.construct_hapID2index_map();
//        this.encoding_haps(); // initialize this.global_haps;
//
//        // TODO: [LEFTOVER]
//        // this.write_global_file_string(
//        //    "/home/lmak/Documents/v0.7_test/Both_V2_Plus_GC_800/gc.inter_freq.txt", false);
//
//    }

    public void construct_hapID2index_map(){
        this.hapID2index = new HashMap<String, Integer>();
        for (int h = 0; h < this.num_global_hap; h++) {
            this.hapID2index.put(this.hap_IDs[h], h);
        }
    }

    /**
     * The HapID is naturally coded by 0 and 1 strings. For very long haplotypes, we hope to convert
     * them in to an integer of base 16. To indicate it is a haplotype, we add an "h" before the 
     * integer 
     */
    public void recode_HapIDs_to_base16() {
    	// check if the current hap-IDs are indeed 0/1 coded
    	for(int h=0;h<this.num_global_hap;h++) {
    		for(int locus=0;locus<this.num_loci;locus++) {
    			if(this.hap_IDs[h].charAt(locus)!='1' && this.hap_IDs[h].charAt(locus)!='0') {
    				System.out.println("Exception: hap-ID "+this.hap_IDs[h]+" is not 1/0 coded."
    						+ "Returned without converting to base 16");
    				return;
    			}
    		}
    	}
    	// after check, convert the IDs:
    	for(int h=0;h<this.num_global_hap;h++) {
    		BigInteger the_integer=new BigInteger(this.hap_IDs[h], 2);
    		this.hap_IDs[h]="h"+the_integer.toString(16);
    	}
    }

    /**
     *  Maps global_haps_string to global_haps (that is coded by integers).
     *
     *  It will be done by invoking the function encoding_alleles in the class LocusAnnotation.
     *  The key algorithm of how to design the encoding will be implemented later by translating 
     *  Michael's MDS code in R to Java TODO [Quan]
     *
     */
    public void encoding_haps(){
        this.global_haps = new double[this.num_global_hap][this.num_loci];
        for (int h = 0; h < this.num_global_hap; h++) {
            for (int l = 0; l < this.num_loci; l++) {
                this.global_haps[h][l] = this.locusInfo[l]
                    .alleles_coding
                    .get(this.global_haps_string[h][l]);

            }
        }
    }


    /**
     *  Prints intermediate haplotype (single-region) guesses to STDOUT.
     */
    public void write_global_stdout() {
        DecimalFormat df = new DecimalFormat("#.####");
        df.setRoundingMode(RoundingMode.CEILING);
        System.out.println("Hap.\tVar. Comp.\tInter. Freq.");
        for (int h = 0; h < this.num_global_hap; h++) {
            System.out.print(this.hap_IDs[h] + "\t");
            for (int l = 0; l < this.num_loci; l++) {
                System.out.print(this.global_haps_string[h][l]);
            }
            System.out.println("\t" + df.format(this.global_haps_freq[h]));
        }
    }


    /**
     *  Output the in global_haplotypes using string alleles.
     *
     *  @param global_hap_output_file
     *  
     */
    public void write_global_file_string(String global_hap_output_file) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(global_hap_output_file));
            bw.write("Hap_ID");
            for (int h = 0; h < this.num_global_hap; h++) {
                bw.write("\t" + this.hap_IDs[h]);
            }
            bw.write("\nFreq");
            for (int h = 0; h < this.num_global_hap; h++) {
                bw.write("\t" + this.global_haps_freq[h]);
            }
            bw.write("\n");
            for (int l = 0; l < this.num_loci; l++) {
                bw.write(this.locusInfo[l].output2string());
                for (int h = 0; h < this.num_global_hap; h++) {
                    bw.write("\t" + this.global_haps_string[h][l]);
                }
                bw.write("\n");
            }
            bw.close();
        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    /**
     *  Output the in global_haplotypes using coded alleles.
     *
     *  During the analysis, e.g., divide-and-conquer, when read by the program, these coded alleles
     *  will become strings and then the next round of coding will be performed.
     *
     *  @param global_hap_output_file
     *  
     */
    public void write_global_file_code(String global_hap_output_file) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(global_hap_output_file));
            bw.write("Hap_ID");
            for (int h = 0; h < this.num_global_hap; h++) {
                bw.write("\t" + this.hap_IDs[h]);
            }

            bw.write("\nFreq");
            for (int h = 0; h < this.num_global_hap; h++) {
                bw.write("\t" + this.global_haps_freq[h]);
            }

            bw.write("\n");
            for (int l = 0; l < this.num_loci; l++) {
                bw.write(this.locusInfo[l].output2string());
                for (int h = 0; h < this.num_global_hap; h++) {
                    // The only difference between this method and "write_global_file_string" is the
                    // line below:
                    bw.write("\t" + this.global_haps[h][l]);
                }
                bw.write("\n");
            }

            bw.close();

        } catch(Exception e) {
            e.printStackTrace();
        }
    }


    /**
     *  Output the in pool frequencies.
     *
     *  @param inpool_hap_output_file
     */
    public void write_inpool(String inpool_hap_output_file) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(inpool_hap_output_file));
            bw.write("Hap_ID");
            for (int h = 0; h < this.num_global_hap; h++) {
                bw.write("\t" + this.hap_IDs[h]);
            }
            bw.write("\n");

            // TODO: [LEFTOVER]
            // System.out.println(this.in_pool_haps_freq.length
            //    + "\t"
            //    + this.in_pool_haps_freq[0].length);

            for (int p = 0; p < this.num_pools; p++) {
                bw.write(this.pool_IDs[p]);
                // TODO (old) [Review]:: Report error! Formerly, this.num_pools.
                for (int h = 0; h < this.num_global_hap; h++) {
                    if (this.in_pool_haps_freq[h].length == 0) {
                        bw.write("\t0");
                    } else {
                        bw.write("\t" + this.in_pool_haps_freq[h][p]);
                    }
                }
                bw.write("\n");
            }
            bw.close();
        } catch(Exception e) {
            e.printStackTrace();
        }
    }


    /**
     *
     *  @param global_hap_output_file
     *  @param in_pool_hap_output_file
     *  @param type
     *  @param append
     */
    public void write2files(
        String global_hap_output_file,
        String in_pool_hap_output_file,
        String type) {

        write_inpool(in_pool_hap_output_file);
        if (type.equals("string")) {
            write_global_file_string(global_hap_output_file);

        } else if(type.equals("code")) {
            write_global_file_code(global_hap_output_file);

        } else {
            System.out.println("Error: the type has to be string or code!");
        }
    }


    /**
     *  Updates the composition of the main global variables when existing haplotypes are removed.
     *
     *  @param list_rem_haps
     *  @param num_rem_haps
     */
    public void remHaps(boolean[] list_rem_haps, int num_rem_haps) {
        int tot_new_haps = this.num_global_hap - num_rem_haps;

        String[][] tmp_global_haps_string = new String[tot_new_haps][this.num_loci];
        String[] tmp_hap_IDs = new String[tot_new_haps];
        double tmp_tot_freq = 0;
        int tmp_count = 0;
        for (int h = 0; h < this.num_global_hap; h++) {
            if (list_rem_haps[h]) {
                continue;
            }
            tmp_global_haps_string[tmp_count] = this.global_haps_string[h].clone();
            tmp_hap_IDs[tmp_count] = this.hap_IDs[h];
            tmp_tot_freq += this.global_haps_freq[h];
            tmp_count++;
        }

        this.global_haps_string = tmp_global_haps_string.clone();
        this.hap_IDs = tmp_hap_IDs.clone();
        double[] tmp_global_haps_freq = new double[tot_new_haps];
        int new_hap_index = 0;
        for (int h = 0; h < this.num_global_hap; h++) {
            if (list_rem_haps[h]) {
                continue;
            }
            tmp_global_haps_freq[new_hap_index] = this.global_haps_freq[h] / tmp_tot_freq;
            new_hap_index++;
        }

        this.global_haps_freq = tmp_global_haps_freq.clone();
        this.num_global_hap = this.global_haps_freq.length;
        this.construct_hapID2index_map();
        this.encoding_haps(); // initialize this.global_haps
    }


    /**
     *  Calculates the average variant frequency, var.-covar. matrix for different primitive loci
     *  Calculates the log-likelihood of this HapConfig object explaining the data.
     *  Updates the global variables mu, sigma, and logL respectively for this HapConfig object.
     */
    public void update_sigma_mu_logL() {
    	
        this.mu = new double[this.num_loci];
        for (int l = 0; l < this.num_loci; l++) {
            for (int h = 0; h < this.num_global_hap; h++) {
                this.mu[l] = this.mu[l] + this.global_haps[h][l] * this.global_haps_freq[h];
            }
        }
        double[][] eta = new double[this.num_loci][this.num_loci];
        for (int l1 = 0; l1 < this.num_loci; l1++) {
             for (int l2 = 0; l2 < this.num_loci; l2++) {
                 for (int h = 0; h < this.num_global_hap; h++) {
                     eta[l1][l2] += (this.global_haps[h][l1]
                        * this.global_haps[h][l2]
                        * this.global_haps_freq[h]);
                  }
             }
        }
        // TODO Update the sigma function to use LDx.
        this.sigma = new double[this.num_loci][this.num_loci];
        for (int q1 = 0; q1 < this.num_loci; q1++) {
             for (int q2 = 0; q2 < this.num_loci; q2++) {
                 this.sigma[q1][q2] = eta[q1][q2] - this.mu[q1] * this.mu[q2];
             }
        }

        this.logL = Algebra.logL_aems(this.sigma, this.mu, this.inpool_site_freqs);
        // Got rid of Algebra.times(this.sigma, this.est_ind_pool),
        // Algebra.times(this.mu, this.est_ind_pool) because everything is in frequencies.
    }


    /**
     *  Returns haplotypes that are above a certain frequency cutoff.
     *
     *  @param freq_cutoff frequency cutoff.
     */
    public HapConfig filter(double freq_cutoff) {
        if (freq_cutoff == 0) {
            return new HapConfig(
                this.global_haps_string,
                this.global_haps_freq,
                this.in_pool_haps_freq,
                this.inpool_site_freqs,
                this.locusInfo, this.num_pools,
                this.hap_IDs,
                this.pool_IDs,
                this.est_ind_pool);
        }

        boolean[] list_rem_haps = new boolean[this.num_global_hap];
        int num_rem_haps = 0;
        for (int h = 0; h < this.num_global_hap; h++) {
            if (this.global_haps_freq[h] < freq_cutoff) {
                list_rem_haps[h] = true;
                num_rem_haps++;
            }
        }
        remHaps(list_rem_haps, num_rem_haps);
        return new HapConfig(
            this.global_haps_string,
            this.global_haps_freq,
            this.in_pool_haps_freq,
            this.inpool_site_freqs,
            this.locusInfo,
            this.num_pools,
            this.hap_IDs,
            this.pool_IDs,
            this.est_ind_pool);

    }
}
