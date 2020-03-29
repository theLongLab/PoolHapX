package PoolHap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;
import java.util.Random;

/**
 * @author Chen Cao 2019-07
 * 
 * 
 */

public class GraphColoring {
    // TODO: [ReconEP]:: refactor variable names to be more sensible (e.g. camelCase for Java).
    // TODO: [Question]:: why are these Vectors instead of ArrayLists in the first place?

    // The index of the genotype when it was added to readinfo_arr_tmp.
    public Vector<Integer> readindex_arr_tmp;

    // The list of all possible genotypes, ArrayList<site=allele; ... site=allele;>
    public Vector<String> readinfo_arr_tmp;
    public HashMap<String, Integer> output_ref_arr ;
    public HashMap<String, String> conf_ref_arr;
    public int num_loci;
    public int num_pools;
    public LocusAnnotation[] locusInfo; // # of loci, note that a locus can be a SNP or a region
    public double[][] inpool_site_freqs; // # of loci x # of pools; Added by Quan Dec. 2018.
    public int num_loci_window;
    public int max_num_gap;
    public String vef_file;
    public double[][] loci_link_count; // (num_loci-1) x 4 (4:	0/0; 0/1/; 1/0; 1/1)
    public double[][] loci_link_freq; // (num_loci-1) x 4 (4:	0/0; 0/1/; 1/0; 1/1) 
    public HashMap<Integer, Integer> loci_index_dict;
    public int bfs_mismatch_tolerance;
    

    // TODO: [LEFTOVER]
    // new GraphColoring(gp.inter_dir + prefix + "_p" + p + ".vef",
    //     gs_var_pos,
    //     gp.inter_dir + prefix + "_p" + p + ".in");

    /**
     *  Graph coloring object constructor with input files. For first round of graph coloring.
     *  TODO: [Javadoc]:: improve description of the constructor for GC round 1.
     *
     *  @param vef variant encoded file file path string.
     *  @param gs_var_pos (required) gold standard variant positions file path string.
     *  @param out_file output file path string.
     *  @throws IOException on input error.
     */
    
    /**
     * @author Chen Cao 2019-08
     * Construct the frequency path for four possible ways: 0/0 0/1 1/0 1/1
     * Using the frequency path to  Evaluate whether the haplotypes ( 2^N )exist 
     * And Estimate the frequency of the haplotypes.
     * The output is used as the initial matrix for AEM 
     * 
     */  
    
    
    
    public GraphColoring(String vef, String gs_var_pos, String out_file, int[][] regions_level_I
    		) throws IOException {
    	this.gc_solver(gs_var_pos, false);
    	this.loci_link_count = new double[this.num_loci-1][4];
    	this.loci_link_freq = new double [this.num_loci-1][4];
    	for (int i = 0; i < this.loci_link_count.length; i++) {
    		for (int j = 0; j < this.loci_link_count[i].length; j++) {
    			this.loci_link_count[i][j] =0.0;
    			this.loci_link_freq[i][j] = 0.0;
    		}
    	}

    	BufferedReader bufferedreader = new BufferedReader(new FileReader(vef));
        String line = "";
        while ((line = bufferedreader.readLine()) != null) {
            line = line.replace("\r", ""); // remove newline characters
//            System.out.println( line);
            String[] line_arr = line.split("\t"); // Read_Name\tPos=Allele;...\t//\tStart\tEnd
            // If the read contains a segregating site (i.e.: has a distinguishing genotype)...
            if (line_arr[1].contains("=")) {
                String str1  = line_arr[1]; // the segregating site information column
                String[] tmp_arr = str1.split(";"); 
                if (tmp_arr.length> 1) {
	                int[] loci_arr = new int[tmp_arr.length];
	                int[] geno_arr = new int[tmp_arr.length] ;
	                for (int i = 0; i < tmp_arr.length; i++) {
	                	String str2  =  tmp_arr[i];
	                	String[] tmp_arr2 = str2.split("=");
	                	loci_arr[i]= Integer.parseInt(tmp_arr2[0] ) ;
	                	geno_arr[i]= Integer.parseInt(tmp_arr2[1] ) ;
	                }
//	                for (int i = 0; i < (geno_arr.length); i++) {
//	                	System.out.println( geno_arr[i]);
//	                }
	                for (int i = 0; i < (geno_arr.length-1); i++) {
//	                	System.out.println( loci_arr[i]);
//	                	System.out.println( loci_arr[i+1]);
	                	if (loci_index_dict.containsKey(loci_arr[i+1]) && 
	                			loci_index_dict.containsKey(loci_arr[i]) ) {
		                	if ((loci_index_dict.get(loci_arr[i+1] ) - loci_index_dict.get(loci_arr[i] )) 
		                		==1 ) {
			                	if ((geno_arr[i]==0)  && (geno_arr[i+1]==0) ) {
			                		this.loci_link_count [this.loci_index_dict.get(loci_arr[i] )][0] += 1.0;
			               		}else if ((geno_arr[i]==0)  && (geno_arr[i+1]==1) ) {
			               			this.loci_link_count [this.loci_index_dict.get(loci_arr[i] )][1] += 1.0;
			               		}else if ((geno_arr[i]==1)  && (geno_arr[i+1]==0) ) {
			               			this.loci_link_count [this.loci_index_dict.get(loci_arr[i] )][2] += 1.0;
			               		}else if ((geno_arr[i]==1)  && (geno_arr[i+1]==1) ) {
			               			this.loci_link_count [this.loci_index_dict.get(loci_arr[i] )][3] += 1.0;
			               		}
		                	}
	                	}
	                }
                }
            }
        }
        bufferedreader.close();
        
    	double lowest_freq_cutoff= 0.03;
    	double loci_min_count = 25.0;
    	
    	for (int i = 0; i < this.loci_link_count.length; i++) {
    		double total_count =0.0;
    		for (int j = 0; j < this.loci_link_count[i].length; j++) {
    			total_count +=   this.loci_link_count[i][j];
    		}
    		if (total_count > loci_min_count) {
    			for (int j = 0; j < this.loci_link_count[i].length; j++) {
    				if (( this.loci_link_count[i][j]/ total_count) >  lowest_freq_cutoff) {
    					this.loci_link_freq[i][j]= this.loci_link_count[i][j]/ total_count;
    				}else {
    					this.loci_link_freq[i][j] = 0.0;
    				}
    			}
    		}else {
    			this.loci_link_freq[i][0] = this.loci_link_freq[i][1]=
    					this.loci_link_freq[i][2]=this.loci_link_freq[i][3]=0.25;
    		}
    	}
//    	for (int i = 0; i < this.loci_link_count.length; i++) {
//    		System.out.print(this.loci_link_freq[i][0]
//                    + "\t"
//                    + this.loci_link_freq[i][1]
//                    + "\t"
//                    + this.loci_link_freq[i][2]
//                    + "\t"
//    				+ this.loci_link_freq[i][3]
//                    + "\n"
//    				);
//    	}
    }
    
    public GraphColoring(String vef, String gs_var_pos, String out_file) throws Exception {
    	HashMap<Integer, Integer> pos_dict = new HashMap<Integer, Integer>();
        Integer pos_index = 0;

        BufferedReader br = new BufferedReader(new FileReader(gs_var_pos));
        String currLine = br.readLine(); // skip header
        currLine = br.readLine();

        while (currLine != null) {
            pos_dict.put(Integer.parseInt(currLine.split(";")[1]), pos_index);


            pos_index++; // move to next variant position index
            currLine = br.readLine(); // read next line
        }
        br.close();
//        this.loci_index_dict = pos_dict;
        this.num_loci = pos_dict.size(); 
        
        
        BufferedWriter bw = new BufferedWriter(new FileWriter(out_file));
        String tmp ="1";
	    for (int i=0;i< (this.num_loci-1);i++) {
	       		 tmp=tmp+"-1";
	    }
	    tmp=tmp+"\t1";
	    bw.write(tmp+"\n");
        bw.close();
       	 
        
        
    }
    
    
    
    public GraphColoring(String vef, String gs_var_pos, String out_file, int num_pos_window, 
    		int num_gap_window) throws IOException {
        /*
         *  Initialize read indices, read information, genotype dictionary, and reader variables.
         */
        // TODO: [Question]:: why are the variables named "arr" when they're Vector or HashMap
        // objects?
        this.readindex_arr_tmp = new Vector<Integer>();
        this.readinfo_arr_tmp = new Vector<String>();
        this.num_loci_window = num_pos_window;
        this.max_num_gap = num_gap_window;
        this.vef_file = vef;
        
        int count = 0;

        // HashMap<pos=allele;..., count>
        HashMap<String, Integer> geno_dict = new HashMap<String, Integer>();

        // I/O for VEF file.
        BufferedReader bufferedreader = new BufferedReader(new FileReader(vef));
        String line = "";

        // The maximum number of times a genotype can be counted in a single pool.
        // TODO: (old) [Question]::  Why does this exist?
        int max_num_geno = 30000;


        /*
         *  Read through VEF file.
         */
        while ((line = bufferedreader.readLine()) != null) {
            line = line.replace("\r", ""); // remove newline characters
            String[] line_arr = line.split("\t"); // Read_Name\tPos=Allele;...\t//\tStart\tEnd

            // TODO: [ReconEP]:: extract the if-else logic below into a separate method.
            // If the read contains a segregating site (i.e.: has a distinguishing genotype)...
            if (line_arr[1].contains("=")) {
                String tmp_geno = line_arr[1]; // the segregating site information column

                // If this combination of alleles hasn't been recorded yet, add to dictionary.
                if (!geno_dict.containsKey(tmp_geno)) {
                    geno_dict.put(tmp_geno, 1);

                    // The index of readinfo_arr_tmp that corresponds to this genotype.
                    this.readindex_arr_tmp.add(count);
                    count = count + 1; // TODO: (minor) [Question]:: change to += for consistency?
                    this.readinfo_arr_tmp.add(line_arr[1]);

                // Else, it has already been recorded, remove and re-add with count + 1.
                } else {
                    int tmp_num = geno_dict.get(tmp_geno); // current count of site combination

                    // TODO: [Question]:: do we need to remove and re-add? Wouldn't the following
                    // auto-update and eliminate most of the if-else outside of guarding against the
                    // maximum?
                    // map.put(key, map.getOrDefault(key, 0) + 1)
                    geno_dict.remove(tmp_geno);
                    geno_dict.put(tmp_geno, (tmp_num + 1));

                    // If site combination count is less than the stated max, re-add to dictionary.
                    if (geno_dict.get(tmp_geno) < max_num_geno) {
                        this.readindex_arr_tmp.add(count);
                        count = count + 1;
                        this.readinfo_arr_tmp.add(line_arr[1]);
                    }
                }
            }
        }
        bufferedreader.close();
        /*
         *  Solve and produce haplotype configurations.
         */
        this.gc_solver(gs_var_pos, true);
        this.fileOut(out_file);
    }
    
	public static String link_str(String x, String y,ArrayList<Integer >  pos_arr  )
			throws IOException {
		String[] x_char = new String [x.length()];
		String[] y_char = new String [y.length()];
		for (int i=0; i< x_char.length;i++) {
			x_char[i]= x.substring(i,i+1);
		}
		for (int i=0; i< y_char.length;i++) {
			y_char[i]= y.substring(i,i+1);
		}
		for (int i=0; i< x_char.length;i++) {
			for (int j=0; j< pos_arr.size();j++) {
				if (i== pos_arr.get(j)) {
					x_char[i] = y_char[j];
				}
			}
		}
		String tmp ="";
		for (int i=0; i< x_char.length;i++) {
			tmp=tmp+ x_char[i];
		}
		return tmp ;
	}
    /**
     * @author Chen Cao 2019-09
     * For the adjacent level I regions, mismatches tolerance (num of 
     * mismatches bases ) is constructed regions linked in the process
     * of breadth first search.
     * 
     */  	
    
    public String [] strmatch (String x, String y, String z, 
    		int num, int num_mismatch_cutoff)  throws IOException {
		
    	ArrayList<Integer >  pos_arr = new ArrayList<Integer>();

    	String overlap_x= x.substring(x.length()-num);
    	String overlap_xz= z.substring(0, num);
    	int distance=0; 
    	
    	for (int i = 0; i < overlap_x.length(); i++) {
			if (overlap_x.charAt(i) != overlap_xz.charAt(i)) {
				distance++;
				pos_arr.add(x.length()- num+i);
			}
		}
    	
    	String overlap_y= y.substring(0, z.length()-num);
    	String overlap_yz= z.substring( num);
    	for (int i = 0; i < overlap_y.length(); i++) {
			if (overlap_y.charAt(i) != overlap_yz.charAt(i)) {
				distance++;
				pos_arr.add(x.length()+ i);
			}
		}
    	if (distance > num_mismatch_cutoff) {
    		String[] com_str2 = new String [ 0];
    		return com_str2;
    	}else {
    		int haps_2n = (int) Math.pow(2, distance); 
        	String[] com_str = new String [ haps_2n];
        	for (int h = 0; h < haps_2n; h++) {
        		com_str[h]= x+y;
        	}
        	String[] sub_str = new String [ haps_2n];
        	for (int h = 0; h < haps_2n; h++) {
        		String curr_ID = "";
            	String vc_str = Integer.toBinaryString(h);
            	int l = vc_str.length();
            	for (int locus = 0; locus < distance - l ; locus++) {
            		vc_str = "0"+ vc_str;
                }
            	sub_str[h] = vc_str;
        	}
    		for (int h = 0; h < haps_2n; h++) {
    			com_str[h]= link_str(com_str[h], sub_str[h], pos_arr );
    		}
    		return com_str;
    	}
    	
    }
    
    
    public boolean reg_strmatch (String x, String y, String z, 
    		int num, int num_mismatch_cutoff)  throws IOException {
		
    	ArrayList<Integer >  pos_arr = new ArrayList<Integer>();

    	String overlap_x= x.substring(x.length()-num);
    	String overlap_xz= z.substring(0, num);
    	int distance=0; 
    	for (int i = 0; i < overlap_x.length(); i++) {
			if (overlap_x.charAt(i) != overlap_xz.charAt(i)) {
				distance++;
				pos_arr.add(x.length()- num+i);
			}
		}
    	String overlap_y= y.substring(0, z.length()-num);
    	String overlap_yz= z.substring( num);
    	for (int i = 0; i < overlap_y.length(); i++) {
			if (overlap_y.charAt(i) != overlap_yz.charAt(i)) {
				distance++;
				pos_arr.add(x.length()+ i);
			}
		}
    	if (distance > num_mismatch_cutoff) {
    		return false;
    	}else {
    		return true;
    	}
    }
    
    public boolean strmatch (String x, String y,  int num) throws IOException {
    	
    	if (num <1) {
    		return true;
    	}
    	String overlap_x= x.substring(x.length()-num);
    	String overlap_y= y.substring(0, num);
    	if (overlap_x.equals(overlap_y)) {
    		return true;
    	}else {
    		return false;
    	}
    }
    
    
    public String  strcombine (String x, String y, int num) throws IOException {
    	if (num==0) {
    		return x+y;
    	} else {
    		return x+y.substring(num);
    	}
    }
    
    public double  minvalue (double x, double y) throws IOException {
    	if (x>= y) {
    		return y;
    	}else {
    		return x;
    	}
    }
    
    public GraphColoring(int start, int end,  HapConfig[] level_1,  HapConfig[] level_2,  String gs_var_pos, 
    		int [][] level_1_region , int [][] level_2_region , int mismatch_tolerance, int reg_ini_region,
    		HashMap<Integer, String> index_var_prefix_dict , String prefix,int num_regeions_merge ) 
    				throws IOException {
    	this.bfs_mismatch_tolerance= mismatch_tolerance;
    	int level_1_first =-1;
    	int level_1_second =-1;
    	for (int i=0; i< level_1_region.length; i++ ) {
    		if (level_1_region[i][0]== start) {
    			level_1_first = i;
    		}
    		if (level_1_region[i][1]== end) {
    			level_1_second = i;
    		}
    	}
    	
    	if (( level_1_first == -1) ||  ( level_1_second == -1)  ){
    		System.out.println("Can not detect the correct Level III regions to generate the initialized "
    				+ "haplotypes for L0L1 regression.");
    		System.exit(0);
    	}	
    	
    	ArrayList<String >  bfs_haps= new ArrayList<String>();
    	ArrayList<Integer >  haps_end_pos= new ArrayList<Integer>();
    	int start_index =0;
    	int end_index=0;
    	for (int i = 0; i  < level_1[level_1_first].global_haps_string.length; i++) {
    		String tmp = "";
    		for (int j = 0; j  < level_1[level_1_first].global_haps_string[i].length; j++) {
    			tmp=tmp+ level_1[level_1_first].global_haps_string[i][j];
    		}
    		bfs_haps.add(tmp);
    		haps_end_pos.add(level_1_region[level_1_first][1]);
    		end_index++;
    	}
    	
    	
    	for (int i = level_1_first ; i  <  level_1_second; i++) {   
    		HashSet<String > hap_dict = new HashSet <String > ();
    		
    		
    		ArrayList<String  >  seg_haps= new ArrayList<String >();
    		int end_pos= haps_end_pos.get(end_index-1);
    		int start_pos= level_1_region[i+1][0];
    		int fix_end_index= end_index;
    		
    		for (int j = 0; j  < level_1[i+1].global_haps_string.length; j++) {
    			String tmp = "";
    			for (int k = 0; k  < level_1[i+1].global_haps_string[j].length; k++) {
    				tmp=tmp+ level_1[i+1].global_haps_string[j][k];
    			}
    			seg_haps.add(tmp );
    		}
    		
    		ArrayList<String  >  stick_haps= new ArrayList<String >();
    		for (int j = 0; j  < level_2[i].global_haps_string.length; j++) {
    			String tmp = "";
    			start_pos= level_2_region[i][0];
    			for (int k = 0; k  < level_2[i].global_haps_string[j].length; k++) {
    				tmp=tmp+ level_2[i].global_haps_string[j][k];
    			}
    			stick_haps.add(tmp );
    		}
    		
    		for (int j = start_index; j  < fix_end_index; j++) {
    			for (int k = 0; k  < seg_haps.size(); k++) {
    				for (int l = 0; l  < stick_haps.size(); l++) {
    					if ( reg_strmatch(bfs_haps.get(j), seg_haps.get(k) ,stick_haps.get(l), 
	    						end_pos- start_pos+1, this.bfs_mismatch_tolerance)) {
    						String tmp_hap =   bfs_haps.get(j)+seg_haps.get(k);
    						if (!hap_dict.contains(tmp_hap)) {
    							hap_dict.add(tmp_hap);
			    				bfs_haps.add( tmp_hap);
			    				haps_end_pos.add(level_1_region[i+1][1]  ); 
			    				end_index++;
			    				break;
    						}
//    						tmp_hap = tmp_hap.substring(0,start_pos)+ stick_haps.get(l)
//    							+tmp_hap.substring(start_pos  + stick_haps.get(l).length());
//    						if (!hap_dict.contains(tmp_hap)) {
//    							hap_dict.add(tmp_hap);
//			    				bfs_haps.add( tmp_hap);
//			    				haps_end_pos.add(level_1_region[i+1][1]  ); 
//			    				end_index++;
//    						}
    						
	    				}
    				}
    			}
    		}
    		start_index =fix_end_index;
    		seg_haps.clear();
    		stick_haps.clear();
    	}
    	
    	
    	ArrayList<String >  ini_haps= new ArrayList<String>();
    	
    	for (int i = start_index; i  < end_index; i++) {
        	ini_haps.add(bfs_haps.get(i)); 
        }

    	String outfile = prefix+ "/regression_level_1_region_"+ Integer.toString(reg_ini_region)+
    			".potential.haps";
    	
    	BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
        bw.write("Hap_ID");
        for (int h = 0; h < ini_haps.size(); h++) {
            bw.write("\th" + Integer.toString(h));
        }

        bw.write("\nFreq");
        for (int h = 0; h < ini_haps.size(); h++) {
            bw.write("\t" + 1/ ((double )ini_haps.size() ) );
        }
        bw.write("\n");
        for (int l = 0; l < end-start+1; l++) {
            bw.write(index_var_prefix_dict.get(start+l) );
            for (int h = 0; h < ini_haps.size(); h++) {
                bw.write("\t" + ini_haps.get(h).substring(l,l+1));
            }
            bw.write("\n");
        }
        bw.close();
    }
    
    /**
     *  Chen: apply Breadth-First-Search (incluing Pruning )to replace gc.
     *  TODO: [Javadoc]:: .
     *
     *  @param level_1 all level 1 haplotype configurations
     *  @param level_2 all level 2 haplotype configurations
     *  @param gs_var_pos (required) gold standard variant positions file path string
     *  @param virtual_cov_link_gc // Do not need this parameter any more. 
     *  @throws IOException on input error.
     *  Adjacent level I regions are linked using the info comes from
     *  level II.
     *  For the adjacent level I regions, mismatches tolerance (num of 
     * 	mismatches bases ) is constructed regions linked in the process
     * 	of breadth first search.
     */
    
// Regression Initialization     
    public GraphColoring(int start, int end,  HapConfig[] level_1,  HapConfig[] level_2,  String gs_var_pos, 
    		int [][] level_1_region , int [][] level_2_region , int mismatch_tolerance, int reg_ini_region,
    		HashMap<Integer, String> index_var_prefix_dict , String prefix) throws IOException {
    	this.bfs_mismatch_tolerance= mismatch_tolerance;
    	int level_1_first =-1;
    	int level_1_second =-1;
    	int level_2_paste =-1;
    	for (int i=0; i< level_1_region.length; i++ ) {
    		if (level_1_region[i][0]== start) {
    			level_1_first = i;
    		}
    		if (level_1_region[i][1]== end) {
    			level_1_second = i;
    		}
    	}
    	for (int i=0; i< level_2_region.length; i++ ) {
    		if ((level_2_region[i][0]<= level_1_region[level_1_first][1] ) 
    				&& (level_2_region[i][1]>= level_1_region[level_1_second][0] )) {
    			level_2_paste =i;
    		}
    	}
    	if (( level_1_first == -1) ||  ( level_1_second == -1) || ( level_2_paste==-1 ) ){
    		System.out.println("Can not detect the correct Level III regions to generate the initialized "
    				+ "haplotypes for L0L1 regression.");
    		System.exit(0);
    	}	
    	
    	ArrayList<String >  ini_haps= new ArrayList<String>();
    	if ( level_1_first == level_1_second) {
    		
    		String [] first_haps = new String [level_1[level_1_first].global_haps_string.length];
    		for (int i=0; i<first_haps.length;i++ ) {
    			String tmp = "";
    			for (int j =0; j <level_1[level_1_first].global_haps_string[i].length;j++ ) {
    				tmp=tmp+ level_1[level_1_first].global_haps_string[i][j]; 
    			}
    			ini_haps.add(tmp);
    		}
    		
    	}else {
//    		System.out.println("len:\t"+ level_1[level_1_first].global_haps_string.length);
    		String [] first_haps = new String [level_1[level_1_first].global_haps_string.length];
    		String [] second_haps = new String [level_1[level_1_second].global_haps_string.length];
    		String [] paste_haps = new String [level_2[level_2_paste].global_haps_string.length];
    		for (int i=0; i<first_haps.length;i++ ) {
    			String tmp = "";
    			for (int j =0; j <level_1[level_1_first].global_haps_string[i].length;j++ ) {
    				tmp=tmp+ level_1[level_1_first].global_haps_string[i][j]; 
    			}
    			
    			first_haps[i] = tmp;
    			
    		}
    		
    		for (int i=0; i<second_haps.length;i++ ) {
    			String tmp = "";
    			for (int j =0; j <level_1[level_1_second].global_haps_string[i].length;j++ ) {
    				tmp=tmp+ level_1[level_1_second].global_haps_string[i][j]; 
    			}
    			second_haps[i] = tmp;
//    			System.out.println("sec\t"+ tmp);
    		}
    		
    		for (int i=0; i<paste_haps.length;i++ ) {
    			String tmp = "";
    			for (int j =0; j <level_2[level_2_paste].global_haps_string[i].length;j++ ) {
    				tmp=tmp+ level_2[level_2_paste].global_haps_string[i][j]; 
    			}
    			paste_haps[i] = tmp;
//    			System.out.println("paste\t"+ tmp);
    		}
    		int overlap = level_1_region[level_1_first][1] -level_2_region[level_2_paste][0] +1;
    		for (int i=0; i<first_haps.length;i++ ) {
    			for (int j =0; j <second_haps.length;j ++ ) {
    				boolean ok_paste =false;
    				for (int k =0; k <paste_haps.length;k++ ) {
//    					System.out.println("?????\t"+ first_haps[i]+"\t"+overlap );
    					if (reg_strmatch( first_haps[i], second_haps[j], paste_haps[k], 
    							overlap,this.bfs_mismatch_tolerance )) {
    						ok_paste =true;
    						break;
    					}
    				}
    				if (ok_paste) {
    					ini_haps.add(first_haps[i]+ second_haps[j]);
    				}
    			}
    		}
    	}
    	
    	
    	

    	String outfile = prefix+ "/regression_level_1_region_"+ Integer.toString(reg_ini_region)+
    			".potential.haps";
    	
    	BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
        bw.write("Hap_ID");
        for (int h = 0; h < ini_haps.size(); h++) {
            bw.write("\th" + Integer.toString(h));
        }

        bw.write("\nFreq");
        for (int h = 0; h < ini_haps.size(); h++) {
            bw.write("\t" + 1/ ((double )ini_haps.size() ) );
        }
        bw.write("\n");
        for (int l = 0; l < end-start+1; l++) {
            bw.write(index_var_prefix_dict.get(start+l) );
            for (int h = 0; h < ini_haps.size(); h++) {
                bw.write("\t" + ini_haps.get(h).substring(l,l+1));
            }
            bw.write("\n");
        }
        
        bw.close();
    }
    
    public GraphColoring( int curr_start, int curr_end, int [][] upper_level,   HapConfig[] level_2, 
    		String gs_var_pos, int [][] level_2_region , int mismatch_tolerance, int reg_region,
    		int regression_level,  HashMap<Integer, String> index_var_prefix_dict , String prefix, 
    		int num_regeions_merge, int number_maximum_selected_haplotypes) 
    				throws IOException {
    	int upper_level_first =-1;
    	int upper_level_second =-1;
    	
    	
    	for (int i=0; i< upper_level.length; i++ ) {
    		if (upper_level[i][0]== curr_start ) {
    			upper_level_first = i;
    		}
    		if (upper_level[i][1]== curr_end ) {
    			upper_level_second = i;
    		}
    	}
    	int [] level_2_paste = new int [upper_level_second- upper_level_first ];
    	
    	for (int u=upper_level_first; u< upper_level_second;u++ ) {
	    	for (int i=0; i< level_2_region.length; i++ ) {
	    		
	    		if ((level_2_region[i][0]<= upper_level[u][1] ) 
	    				&& (level_2_region[i][1]>= upper_level[u+1][0] )) {
	    			level_2_paste[u-upper_level_first ] =i;
	    		}
	    	}
    	}
    	

    	if (( upper_level_first == -1) ||  ( upper_level_second == -1) ){
    		System.out.println("Can not detect the correct level "+  Integer.toString(regression_level)
    				+ " regression region: "  + 
    				Integer.toString(reg_region)+"\n");
    		System.exit(0);
    	}	
    	
    	
    	ArrayList<String >  upper_level_haps= new ArrayList<String>();
    	ArrayList<Double >  upper_level_haps_freq= new ArrayList<Double >();
    	ArrayList<String >  combine_haps= new ArrayList<String>();
    	ArrayList<String >  tmp_combine_haps= new ArrayList<String>();
    	combine_haps.add("");
    	for (int u= upper_level_first; u<=  upper_level_second;u++ ) {
    		upper_level_haps.clear();
    		upper_level_haps_freq.clear();
    		tmp_combine_haps.clear();
        	String upper_level_first_file = prefix+ "/regression_level_"+ Integer.toString(regression_level-1)+
        			"_region_"+Integer.toString(u+1)+".regression_out";
        	BufferedReader bufferedreader1 = new BufferedReader(new FileReader(upper_level_first_file));
        	String line = "";
            while ((line = bufferedreader1.readLine()) != null) {
            	if (!line.substring(0, 1).equals("#")) {
	            	line= line.replace("\n", "").replace("\r", "");
	//            	System.out.println(line);
	            	String[] tmp = line.split("\t");
	            	if (Double.parseDouble(tmp[1]) >0 ) {
	            		upper_level_haps.add(tmp[2]); 
	            		upper_level_haps_freq.add(Double.parseDouble(tmp[1])); 
	            	}
            	}
            }
            bufferedreader1.close();
            double [] tmp_freq =new double [upper_level_haps_freq.size()];
            for (int i=0 ;i < tmp_freq.length; i++) {
            	tmp_freq[i]= upper_level_haps_freq.get(i);
            }
            for (int i=0 ;i < tmp_freq.length; i++) {
            	for (int j=i ;j < tmp_freq.length; j++) {
            		if (tmp_freq[i]< tmp_freq[j]) {
            			double tmp_double= tmp_freq[i];
            			tmp_freq[i] = tmp_freq[j];
            			tmp_freq[j]= tmp_double;
            		}
            	}
            }
            double freq_cutoff= 0.0;
            if (tmp_freq.length > number_maximum_selected_haplotypes) {
            	freq_cutoff= tmp_freq[number_maximum_selected_haplotypes-1];
            }
//            System.out.println(freq_cutoff);
            for (int i=0;i< combine_haps.size();i++) {
            	tmp_combine_haps.add(combine_haps.get(i))  ;
            }
            combine_haps.clear();
            for (int i=0; i< tmp_combine_haps.size();i++) {
            	for (int j =0; j < upper_level_haps.size();j++) {
            		if (upper_level_haps_freq.get(j) >= freq_cutoff) {
            			boolean flag= false;
            			if (u== upper_level_first ) {
            				flag =true;
            			}else {
            				int l2_index = level_2_paste[u-upper_level_first-1 ];
            				for (int l2=0; l2< level_2[l2_index].global_haps_string.length;l2++) {
            					String tmp_str="";
            					for (int k=0;k<level_2[l2_index].global_haps_string[l2].length;k++) {
            						tmp_str = tmp_str+ level_2[l2_index].global_haps_string[l2][k];
            					}
            					
            					int overlap = upper_level[u-1][1] -level_2_region[l2_index][0] +1;
//            					System.out.println(upper_level[u][0]+"\t"+upper_level[u][1]+"\t"+overlap);
            					if (reg_strmatch(tmp_combine_haps.get(i),  upper_level_haps.get(j), 
            							tmp_str,overlap, mismatch_tolerance) ) {
            						flag=true;
            						break;
            					}
            					
            				}
            			}
            			if (flag){
            				combine_haps.add(tmp_combine_haps.get(i)+upper_level_haps.get(j) );
//            			System.out.println(tmp_combine_haps.get(i)+upper_level_haps.get(j));
            			}
            		}
            	}
            }
            
    	}
//    	System.out.println(combine_haps.size() );

    	String outfile = prefix+ "/regression_level_"+ Integer.toString(regression_level)+
    			"_region_"+ Integer.toString(reg_region)+ ".potential.haps";
    	
    	BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
        bw.write("Hap_ID");
        for (int h = 0; h < combine_haps.size(); h++) {
            bw.write("\th" + Integer.toString(h));
        }

        bw.write("\nFreq");
        for (int h = 0; h < combine_haps.size(); h++) {
            bw.write("\t" + 1/ ((double )combine_haps.size() ) );
        }
        bw.write("\n");
        for (int l = 0; l < curr_end-curr_start+1; l++) {
            bw.write(index_var_prefix_dict.get(curr_start+l) );
            for (int h = 0; h < combine_haps.size(); h++) {
                bw.write("\t" + combine_haps.get(h).substring(l,l+1));
            }
            bw.write("\n");
        }
        bw.close();
    }
    
    
    
    
    public GraphColoring(int curr_start, int curr_end, int [][] upper_level,   HapConfig[] level_2, 
    		String gs_var_pos, int [][] level_2_region , int mismatch_tolerance, int reg_region,
    		int regression_level,  HashMap<Integer, String> index_var_prefix_dict , String prefix) 
    				throws IOException {
    	int upper_level_first =-1;
    	int upper_level_second =-1;
    	int level_2_paste =-1;
    	for (int i=0; i< upper_level.length; i++ ) {
    		if (upper_level[i][0]== curr_start ) {
    			upper_level_first = i;
    		}
    		if (upper_level[i][1]== curr_end ) {
    			upper_level_second = i;
    		}
    	}
    	
    	for (int i=0; i< level_2_region.length; i++ ) {
    		if ((level_2_region[i][0]<= upper_level[upper_level_first][1] ) 
    				&& (level_2_region[i][1]>= upper_level[upper_level_second][0] )) {
    			level_2_paste =i;
    		}
    	}
    	

    	if (( upper_level_first == -1) ||  ( upper_level_second == -1) || ( level_2_paste==-1 ) ){
    		System.out.println("Can not detect the correct level "+  Integer.toString(regression_level)
    				+ " regression region: "  + 
    				Integer.toString(reg_region)+"\n");
    		System.exit(0);
    	}	
    	
    	
    	ArrayList<String >  final_haps= new ArrayList<String>();
    	ArrayList<String >  upper_level_first_haps= new ArrayList<String>();
    	ArrayList<String > upper_level_second_haps= new ArrayList<String>();
    	String upper_level_first_file = prefix+ "/regression_level_"+ Integer.toString(regression_level-1)+
    			"_region_"+Integer.toString(upper_level_first+1)+".regression_out";
    	String upper_level_second_file = prefix+ "/regression_level_"+ Integer.toString(regression_level-1)+
    			"_region_"+Integer.toString(upper_level_second+1)+".regression_out";
    	BufferedReader bufferedreader1 = new BufferedReader(new FileReader(upper_level_first_file));
    	String line = "";
        while ((line = bufferedreader1.readLine()) != null) {
        	line= line.replace("\n", "").replace("\r", "");
        	String[] tmp = line.split("\t");
        	if (Double.parseDouble(tmp[1]) >0 ) {
        		upper_level_first_haps.add(tmp[2]); 
        	}
        }
        bufferedreader1.close();
        BufferedReader bufferedreader2 = new BufferedReader(new FileReader(upper_level_second_file));
        while ((line = bufferedreader2.readLine()) != null) {
        	line= line.replace("\n", "").replace("\r", "");
        	String[] tmp = line.split("\t");
        	if (Double.parseDouble(tmp[1]) >0 ) {
        		upper_level_second_haps.add(tmp[2]); 
        	}
        }
        bufferedreader2.close();
    	if ( upper_level_first == upper_level_second) {	
    		for (int i=0; i<upper_level_first_haps.size();i++ ) {	
    			final_haps.add(upper_level_first_haps.get(i));
    		}
    		
    	}else {
//    		System.out.println("len:\t"+ level_1[level_1_first].global_haps_string.length);
    		String [] paste_haps = new String [level_2[level_2_paste].global_haps_string.length];
    		for (int i=0; i<paste_haps.length;i++ ) {
    			String tmp = "";
    			for (int j =0; j <level_2[level_2_paste].global_haps_string[i].length;j++ ) {
    				tmp=tmp+ level_2[level_2_paste].global_haps_string[i][j]; 
    			}
    			paste_haps[i] = tmp;
//    			System.out.println("paste\t"+ tmp);
    		}
    		int overlap = upper_level[upper_level_first][1] -level_2_region[level_2_paste][0] +1;
    		for (int i=0; i<upper_level_first_haps.size();i++ ) {
    			for (int j =0; j <upper_level_second_haps.size();j ++ ) {
    				boolean ok_paste =false;
    				for (int k =0; k <paste_haps.length;k++ ) {
//    					System.out.println("?????\t"+ first_haps[i]+"\t"+overlap );
    					if (reg_strmatch( upper_level_first_haps.get(i), upper_level_second_haps.get(j),
    							paste_haps[k], overlap, mismatch_tolerance )) {
    						ok_paste =true;
    						break;
    					}
    				}
    				if (ok_paste) {
    					final_haps.add(upper_level_first_haps.get(i)+ 
    							upper_level_second_haps.get(j));
    				}
    			}
    		}
    	}
    	
    	
    	

    	String outfile = prefix+ "/regression_level_"+ Integer.toString(regression_level)+
    			"_region_"+ Integer.toString(reg_region)+ ".potential.haps";
    	
    	BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
        bw.write("Hap_ID");
        for (int h = 0; h < final_haps.size(); h++) {
            bw.write("\th" + Integer.toString(h));
        }

        bw.write("\nFreq");
        for (int h = 0; h < final_haps.size(); h++) {
            bw.write("\t" + 1/ ((double )final_haps.size() ) );
        }
        bw.write("\n");
        for (int l = 0; l < curr_end-curr_start+1; l++) {
            bw.write(index_var_prefix_dict.get(curr_start+l) );
            for (int h = 0; h < final_haps.size(); h++) {
                bw.write("\t" + final_haps.get(h).substring(l,l+1));
            }
            bw.write("\n");
        }
        
        bw.close();
    }
    
    
    
    

    public GraphColoring(HapConfig[] level_1,  HapConfig[] level_2,  String gs_var_pos, 
    		int [][] level_1_region , int [][] level_2_region , int mismatch_tolerance
            ) throws IOException {
    	this.bfs_mismatch_tolerance= mismatch_tolerance;
    	ArrayList<String >  bfs_haps= new ArrayList<String>();
    	ArrayList<Double >  haps_min_freq= new ArrayList<Double>();
    	ArrayList<Integer >  haps_end_pos= new ArrayList<Integer>();
    	ArrayList<Integer >  haps_score= new ArrayList<Integer>();
    	int start_index =0;
    	int end_index=0;
    	double aver_freq= 1/ (double) level_1[0].global_haps_freq.length;
    	for (int i = 0; i  < level_1[0].global_haps_string.length; i++) {
    		String tmp = "";
    		for (int j = 0; j  < level_1[0].global_haps_string[i].length; j++) {
    			tmp=tmp+ level_1[0].global_haps_string[i][j];
    		}
    		bfs_haps.add(tmp);
    		haps_min_freq.add(level_1[0].global_haps_freq[i]);
    		haps_end_pos.add(level_1_region[0][1]);
    		if (level_1[0].global_haps_freq[i]>= aver_freq) {
    			haps_score.add(0);
    		}else {
    			haps_score.add(-1);
    		}
    		end_index++;
    	}
    	
    	boolean  one_one_search= true;
//Chen:     Breadth-First-Search
    	//Search Strategy: Level I 1 -> Level I 2 -> Level I 3 -> Level I 4...
    	if (one_one_search) {
	    	for (int i = 0; i  < (level_1_region.length -1); i++) {   
//	    		System.out.println( i );
//	    		System.out.println( end_index-start_index  ); 
	    		ArrayList<String  >  seg_haps= new ArrayList<String >();
	    		ArrayList<Double  >  seg_freq= new ArrayList<Double >();
	    		int end_pos= haps_end_pos.get(end_index-1);
	    		int start_pos= level_1_region[i+1][0];
	    		int fix_end_index= end_index;
	    		
	    		for (int j = 0; j  < level_1[i+1].global_haps_string.length; j++) {
	    			String tmp = "";
	    			for (int k = 0; k  < level_1[i+1].global_haps_string[j].length; k++) {
	    				tmp=tmp+ level_1[i+1].global_haps_string[j][k];
	    			}
	    			seg_haps.add(tmp );
	    			seg_freq.add(level_1[i+1].global_haps_freq[j]); 
	    		}
	    		
	    		ArrayList<String  >  stick_haps= new ArrayList<String >();
	    		for (int j = 0; j  < level_2[i].global_haps_string.length; j++) {
	    			String tmp = "";
	    			start_pos= level_2_region[i][0];
	    			for (int k = 0; k  < level_2[i].global_haps_string[j].length; k++) {
	    				tmp=tmp+ level_2[i].global_haps_string[j][k];
	    			}
	    			stick_haps.add(tmp );
	    		}
	    		
	    		for (int j = start_index; j  < fix_end_index; j++) {
	    			for (int k = 0; k  < seg_haps.size(); k++) {
	    				for (int l = 0; l  < stick_haps.size(); l++) {
	    					String [] com_str = strmatch(bfs_haps.get(j), seg_haps.get(k) ,stick_haps.get(l), 
		    						end_pos- start_pos+1, this.bfs_mismatch_tolerance);
		    				if ((com_str.length!= 0)  && (haps_score.get(j)> -2)){
		    					for (int h = 0; h  < com_str.length; h++) {
			    					bfs_haps.add( com_str[h]);
			    					haps_min_freq.add(minvalue(haps_min_freq.get(j), seg_freq.get(k) )  ); 
			    					haps_end_pos.add(level_1_region[i+1][1]  ); 
			    					end_index++;
			    					if (seg_freq.get(k) >= (1/ (double) seg_haps.size())) {
			    						haps_score.add(haps_score.get(j)  );
			    					}else {
			    						haps_score.add(haps_score.get(j)-1);
			    					}
		    					}
		    				}
	    				}
	    			}
	    		}
	    		start_index =fix_end_index;
	    		seg_haps.clear();
	    		seg_freq.clear();
	    		stick_haps.clear();
	    	}
    	}
    	
    	//Search Strategy: Level I 1 -> Level II 1 -> Level I 1 -> Level II 2...
    	if (!one_one_search) {
	    	for (int i = 0; i  < level_2_region.length; i++) {   
	    		ArrayList<String  >  seg_haps= new ArrayList<String >();
	    		ArrayList<Double  >  seg_freq= new ArrayList<Double >();
	    		for (int j = 0; j  < level_2[i].global_haps_string.length; j++) {
	    			String tmp = "";
	    			for (int k = 0; k  < level_2[i].global_haps_string[j].length; k++) {
	    				tmp=tmp+ level_2[i].global_haps_string[j][k];
	    			}
	    			seg_haps.add(tmp );
	    			seg_freq.add(level_2[i].global_haps_freq[j]); 
	    		}
	    		int end_pos= haps_end_pos.get(end_index-1);
	    		int start_pos= level_2_region[i][0];
	    		int fix_end_index= end_index;
	    		for (int j = start_index; j  < fix_end_index; j++) {
	    			for (int k = 0; k  < seg_haps.size(); k++) {
	    				if ((strmatch(bfs_haps.get(j), seg_haps.get(k) ,end_pos- start_pos+1)) 
	        					&& (haps_score.get(j) > (-1*(i+2) ) ) ){
	    					bfs_haps.add( strcombine ( bfs_haps.get(j), seg_haps.get(k) ,end_pos- start_pos+1) );
	    					haps_min_freq.add(minvalue(haps_min_freq.get(j), seg_freq.get(k) )  ); 
	    					haps_end_pos.add(level_2_region[i][1]  ); 
	    					if (seg_freq.get(k) >= (1/ (double) seg_haps.size())) {
	    						haps_score.add(haps_score.get(j)  );
	    					}else {
	    						haps_score.add(haps_score.get(j)-1);
	    					}
	    					end_index++;
	    				}
	    			}
	    		}
	    		
	    		start_index =fix_end_index;
	    		seg_haps.clear();
	    		seg_freq.clear();
	    		for (int j = 0; j  < level_1[i+1].global_haps_string.length; j++) {
	    			String tmp = "";
	    			for (int k = 0; k  < level_1[i+1].global_haps_string[j].length; k++) {
	    				tmp=tmp+ level_1[i+1].global_haps_string[j][k];
	    			}
	    			seg_haps.add(tmp );
	    			seg_freq.add(level_1[i+1].global_haps_freq[j]); 
	    		}
	
	    		end_pos= haps_end_pos.get(end_index-1);
	    		start_pos= level_1_region[i+1][0];
	    		fix_end_index= end_index;
	    		for (int j = start_index; j  < fix_end_index; j++) {
	    			for (int k = 0; k  < seg_haps.size(); k++) {
	    				if ((strmatch(bfs_haps.get(j), seg_haps.get(k) ,end_pos- start_pos+1)) 
	    					&& (haps_score.get(j) >  (-1*(i+2) ) ) ){
	    					bfs_haps.add( strcombine ( bfs_haps.get(j), seg_haps.get(k) ,end_pos- start_pos+1) );
	    					haps_min_freq.add(minvalue(haps_min_freq.get(j), seg_freq.get(k) )  ); 
	    					haps_end_pos.add(level_1_region[i+1][1]  ); 
	    					if (seg_freq.get(k) >= (1/ (double) seg_haps.size())) {
	    						haps_score.add(haps_score.get(j)  );
	    					}else {
	    						haps_score.add(haps_score.get(j)-1);
	    					}
	    					end_index++;
	    				}
	    			}
	    		}
	    		start_index =fix_end_index;
	    	}
    	}
    	

        this.output_ref_arr = new HashMap<String,Integer>();
        ArrayList<String>  hap_list = new ArrayList<String>(); 
        
        for (int i = start_index; i  < end_index; i++) {
        	if ( (i%1)==0) {
        		hap_list.add(bfs_haps.get(i)); 
        	}
        }


        for (int i = 0; i  < hap_list.size(); i++) {
        	this.output_ref_arr.put(hap_list.get(i) , 1); 
        }
//        
          
        this.gc_solver(gs_var_pos, false);
//    	System.exit(0);
    }
    
    
    /**
     *  2nd round of graph coloring for global haplotype through region linking.
     *  TODO: [Javadoc]:: improve description of this 2nd GC round constructor.
     *
     *  @param level_1 all level 1 haplotype configurations
     *  @param level_2 all level 2 haplotype configurations
     *  @param gs_var_pos (required) gold standard variant positions file path string
     *  @param virtual_cov_link_gc // See property file for the explanation of this parameter
     *  @throws IOException on input error.
     */
    public GraphColoring(HapConfig[] level_1,  HapConfig[] level_2,  String gs_var_pos,
        int virtual_cov_link_gc) throws IOException {
    	this.num_loci_window = Integer.MAX_VALUE;
        this.max_num_gap = 0;
        // Initialize read indices and read information variables.
        this.readindex_arr_tmp = new Vector<Integer>();
        this.readinfo_arr_tmp = new Vector<String>();
        int count = 0;
        // Covert regional haplotypes to VEF format.
        // For each region in all level 1 regions...
        for (int r = 0; r < level_1.length; r++) {
            // For each haplotype in all global haplotypes in region...
            for (int h = 0; h < level_1[r].num_global_hap; h++) {
                String curr_vc = ""; // TODO: [Question]:: what is vc?
                // For each locus in all loci in region
                for (int l = 0; l < level_1[r].num_loci; l++) {
                    curr_vc += (level_1[r].locusInfo[l].start_loc // format read info
                        + "="
                        + level_1[r].global_haps_string[h][l]
                        + ";");
                }
                // TODO: [Question]:: what is ct? Count?
                int hap_ct = (int) (level_1[r].global_haps_freq[h] * virtual_cov_link_gc); // haplotype count(?)
                for (int c = 0; c < hap_ct; c++) { // add reads and corresponding indices to vectors
                    // The index of readinfo_arr_tmp that corresponds to this genotype.
                    this.readindex_arr_tmp.add(count);  // the number of line in the file
                    count++;
                    this.readinfo_arr_tmp.add(curr_vc);
                }
            }
        }
        // For each region in all level 2 regions...
        for (int r = 0; r < level_2.length; r++) {
            // For each haplotype in all global haplotypes in region...
            for (int h = 0; h < level_2[r].num_global_hap; h++) {
                String curr_vc = "";
                // For each locus in all loci in region...
                for (int l = 0; l < level_2[r].num_loci; l++) {
                    curr_vc += (level_2[r].locusInfo[l].start_loc // format read info
                    + "="
                    + level_2[r].global_haps_string[h][l]
                    + ";");
                }
                int hap_ct = (int) (level_2[r].global_haps_freq[h] * virtual_cov_link_gc); // haplotype count(?)
                for (int c = 0; c < hap_ct; c++) { // add reads and corresponding indices to vectors
                    // The index of readinfo_arr_tmp that corresponds to this genotype.
                    this.readindex_arr_tmp.add(count);
                    count++;
                    this.readinfo_arr_tmp.add(curr_vc);
                }
            }
        }
        System.out.println("There are " + count
            + " individual fragments (reads or regional haplotypes) in the dataset.");

        this.gc_solver(gs_var_pos, true);
    }


    /**
     *  TODO: [Question]:: what does this mean?
     *  Source_Path or source_path
     *
     *  Connects variants with graph colouring.
     *  TODO: [Javadoc]:: improve description of solver method.
     *  TODO: [ReconEP]:: this method is way too long, simply not maintainable or testable, refactor
     *  into multiple smaller helper methods.
     *
     *  @param gs_var_pos (required) gold standard variant positions file path string.
     *  @throws IOException on input error.
     */
    
    
    public void gc_solver(String gs_var_pos , boolean run_gc) throws IOException {
        /*
         *  Initialize variables.
         */
        // HashMap<Seg_Site, Index>
        HashMap<Integer, Integer> pos_dict = new HashMap<Integer, Integer>();
        Integer pos_index = 0;

        /*
         *  Read through gold standard variant positions file to load into dictionary and count
         *  total number of loci.
         *
         *  also this needs to be done in SiteInPoolFreqAnno as well, refactor to a static method?
         *  retain number of loci as a global variable?
         */
        // Read into file.
        BufferedReader br = new BufferedReader(new FileReader(gs_var_pos));
        String currLine = br.readLine(); // skip header
        currLine = br.readLine();

        // Read each line into a dictionary.
        while (currLine != null) {
            pos_dict.put(Integer.parseInt(currLine.split(";")[1]), pos_index);

            // TODO: [LEFTOVER]
            // System.out.println(currLine.split(";")[1] + "\t" +  pos_index);

            pos_index++; // move to next variant position index
            currLine = br.readLine(); // read next line
        }
        br.close();
        this.loci_index_dict = pos_dict;
        this.num_loci = pos_dict.size(); // number of loci in gold standard variants file
        

        /*
         *  Read through gold standard variant positions file again to load loci info into a matrix.
         */
        // Read into file.
        br = new BufferedReader(new FileReader(gs_var_pos));
        currLine = br.readLine(); // skip header
        currLine = br.readLine();

        // Initialize variables.
        int loci_index = 0;
        this.num_pools = currLine.split("\t").length - 1;
        this.locusInfo = new LocusAnnotation[this.num_loci];
        this.inpool_site_freqs = new double[this.num_loci][this.num_pools];
        
        while (currLine != null) {
            String[] tmp = currLine.split("\t");

            // The first column is the locus-info.
            this.locusInfo[loci_index] = new LocusAnnotation(tmp[0]);

            // For each pool, load locus frequency into matrix.
            for (int p = 0; p < this.num_pools; p++) {
                this.inpool_site_freqs[loci_index][p] = Double.parseDouble(tmp[p + 1]);
            }

            loci_index++; // move to next locus index
            currLine = br.readLine(); // read next line
        }
        
        
        if (  run_gc ==false) {
        	return;
        }

        // Read through file.
       

        br.close();

        // TODO: [LEFTOVER]
        // System.out.println("There are " + this.locusInfo.length + " positions.");

        /*
         *  Randomize
         *  TODO: [IMPORTANT]:: this entire next section makes no sense to me, it feels like there
         *  are a lot of inefficiencies and overhead that could be avoided by not using Vectors
         *  in the first place? We could literally replace the next 30 lines with like 4 lines of
         *  code if some variables were initialized as ArrayLists.
         */
        // ArrayList<pos=allele;> in random order
        Vector<Integer> readindex_arr = new Vector<Integer>();

        // ArrayList<original index of pos=allele;> TODO: (old) [Review] Confirm!
        Vector<String> readinfo_arr = new Vector<String>();

        // TODO: [Question]:: why not .toArray()?
        // Create new array from read index Vector.
        int[] index_arr_tmp = new int[this.readindex_arr_tmp.size()];
        for (int k = 0; k < index_arr_tmp.length; k++) {
            index_arr_tmp[k] = k;
        }

        // TODO: [Question]:: why not convert from array to arraylist?
        // TODO: [ReconEP]:: rename variables to something that intuitively makes sense...
        // e.g. ArrayList<Integer> list = ArrayList<Integer>(Arrays.asList(index_arr_tmp));
        // Create new ArrayList from read index array.
        ArrayList<Integer> list = new ArrayList<Integer>();
        for (int i = 0; i < index_arr_tmp.length; i++) {
            list.add(index_arr_tmp[i]);
        }

        // A list of the indices corresponding to each read, in random order.
        // Create new array the length of the read index Vector.
        int[] index_arr = new int[this.readindex_arr_tmp.size()];

        // index_arr_tmp, list, and index_arr are copies of this.readindex_arr_tmp so far.
        // Shuffle read index ArrayList.
        Collections.shuffle(list,new Random(19880122));

        // Now, index_arr_tmp, index_arr are copies of this.readindex_arr_tmp. List is now a
        // randomized version of this.readindex_arr_tmp.
        // TODO: [Question]:: why can't we clone the the shuffled list as index_arr? Is iterating
        // through it faster?
        // Honestly
        // e.g. int[] index_arr = list.toArray(new int[list.size()]);
        // Iterate through randomized read index ArrayList to create new read index array that is
        // randomized.
        Iterator<Integer> ite = list.iterator();
        int tmp_i = 0;
        while (ite.hasNext()) {

            // TODO: [LEFTOVER]
            // System.out.println(ite.next().toString()+", ");

            index_arr[tmp_i] = ite.next();
            tmp_i++;
        }

        // Now, index_arr_tmp is a copy of this.readindex_arr_tmp. List, index_arr are now
        // randomized versions of this.readindex_arr_tmp.

        // TODO: [LEFTOVER]
        // for (int i = 0; i < index_arr.length; i++) {
        //     System.out.println(index_arr[i]);
        // }
        // index_arr is a randomized list of indices corresponding to the genotypes from the VEF
        // file.

        // Iteratively add the randomized read index and non-randomized read information arrays
        // to the previously defined new Vectors.
        // TODO: [Question]:: I still don't understand why we need Vectors.
        for (int i = 0; i < this.readindex_arr_tmp.size(); i++) {
            // TODO: [Question]:: why not this.readinfo_arr_tmp?
            readinfo_arr.add(readinfo_arr_tmp.get(index_arr[i]));
            readindex_arr.add(this.readindex_arr_tmp.get(index_arr[i]));
        }
        // readinfo_arr contains the randomized order list of the genotypes from the VEF file.
        // readindex_arr contains the original indices corresponding to the genotypes from the VEF
        // file i.e.: the order they were read in.

        // TODO: [LEFTOVER]
        // for (int i = 0; i < readindex_arr.size(); i++) {
        // 	System.out.println(readinfo_arr.get(i));
        // }


        /*
         *  Set up matrices?
         *  TODO: [Question]:: clarify/understand this section better.
         */
        // Genotype index: ArrayList<ArrayList<Index_SS>>
        ArrayList<ArrayList<Integer>> read_pos_2D_arr= new ArrayList<ArrayList<Integer>>();

        // Genotype allele: ArrayList<ArrayList<Allele_SS>>
        ArrayList<ArrayList<String>> read_geno_2D_arr= new ArrayList<ArrayList<String>>();

        // TODO: [LEFTOVER]
        // System.out.println(readinfo_arr);

        // For each read in the read info array....
        for (int i = 0; i < readinfo_arr.size(); i++) {
            String tmp_str = readinfo_arr.get(i);
            String [] tmp_str_arr= tmp_str.split(";"); // split loci in the read
            ArrayList<String> i_geno_arr = new ArrayList<String>();
            ArrayList<Integer> i_pos_arr = new ArrayList<Integer>();

            // For each locus in the read...
            for (int k = 0; k < tmp_str_arr.length; k++) {
                String [] tmp2_arr = tmp_str_arr [k].split("=");
                int pos = Integer.parseInt(tmp2_arr[0]);
                i_pos_arr.add(pos_dict.get(pos)); // add position index to position row array
                i_geno_arr.add(tmp2_arr[1]); // add genotype to genotype row array
            }

            read_pos_2D_arr.add(i_pos_arr); // add row to read variant positions matrix
            read_geno_2D_arr.add(i_geno_arr); // add row to read genotype matrix

        }

        // TODO: [LEFTOVER]
        // System.out.println(read_pos_2D_arr);
        // System.out.println(read_geno_2D_arr);

        // ArrayList<Genotype, ArrayList<Index_Genotype>>
        ArrayList<ArrayList<Integer>> readadj_arr = new ArrayList<ArrayList<Integer>>();

        // For each possible genotype...
        for (int i = 0; i < read_pos_2D_arr.size(); i++) {

            // TODO: [LEFTOVER]
            // if (i % 500 == 0) {
            //     System.out.println(i + " fragments have been processed.");
            // }

            // TODO: [Question]:: what's this, adjacent?
            ArrayList<Integer> adj_arr = new ArrayList<Integer>();

            // ...comparing it to all other genotypes...
            for (int j = 0; j < read_pos_2D_arr.size(); j++) {
                if (i != j) {
                    boolean IsConnect = false;
                    for (int k = 0; k < read_pos_2D_arr.get(i).size(); k++) {
                        for (int l = 0; l < read_pos_2D_arr.get(j).size(); l++) {
                            // If pos(geno_i, index_k) == pos(geno_j, index_l) and the alleles are
                            // the same, the two reads can connect.
                            if ((read_pos_2D_arr.get(i).get(k) == read_pos_2D_arr.get(j).get(l))
                                && (!read_geno_2D_arr.get(i)
                                    .get(k)
                                    .equals(read_geno_2D_arr.get(j).get(l)))) {

                                IsConnect = true;
                            }
                        }
                    }

                    if (IsConnect){
                        // Add index of geno_j to the list of possible connects for geno_i.
                        adj_arr.add(readindex_arr.get(j));
                    }
                }
            }

            readadj_arr.add(adj_arr);

        }

        // TODO: [LEFTOVER]
        // System.out.println(
        //     "Finished identifying all possible conflicts between genotype fragments.");
        //
        // System.out.println(readadj_arr);

        int max_color = readindex_arr.size(); // only if there are all 2^loci genotypes present
        ArrayList<Integer> read_color_arr = new ArrayList<Integer>();
        ArrayList<HashSet<Integer>> nb_color_arr = new ArrayList<HashSet<Integer>>();

        // TODO: [LEFTOVER]
        // nb_color_arr.add(new HashSet());

        for (int i = 0; i < readadj_arr.size(); i++) {
            nb_color_arr.add(new HashSet<Integer>());
            read_color_arr.add(-1);
        }

        read_color_arr.set(0, 0); // set the first colour as the genotype at index 0

        for (int i = 0; i < readadj_arr.get(0).size(); i++) {
            int index = readadj_arr.get(0).get(i);

            // Add colour 0 as a possible colour to all genotypes that can connect with the genotype
            // index 0.
            nb_color_arr.get(index).add(0);
        }

        int real_max_color = 0;
        ArrayList<HashSet<String>> color_geno_set_arr = new ArrayList<HashSet<String>>();
        color_geno_set_arr.add(new HashSet<String>());

        // This is the full genotype of a colour.
        color_geno_set_arr.get(0).add(readinfo_arr.get(0));

        // TODO: [LEFTOVER]
        // System.out.println(color_geno_set_arr);

        // Make the conflict graph. Basically, assign colours to any genotype fragments that can't
        // go together.
        while (true) {
            int max_nb_color = -1;
            int index = -1;

            for (int i = 0; i < readadj_arr.size(); i++) { // for each genotype...
                // If there hasn't been a colour assigned to that genotype...
                if ((read_color_arr.get(i) == -1)
                    // ...and there are other genotypes that can go with it of the same colour...
                    && (nb_color_arr.get(i).size() > max_nb_color)) {

                    index = i;
                    max_nb_color = nb_color_arr.get(i).size();
                }
            }

            if (index == -1) { // If there are no colours left (?) end this step.
                break;
            }
            int color = -1;

            for (int i = 0; i < max_color; i++) { // for each genotype...
                // If the possible colours list for that genotype doesn't contain colour i...
                if (!nb_color_arr.get(index).contains(i)) {

                    // If colour i is smaller than the number of available colours.
                    if (i < color_geno_set_arr.size()) {
                        if (!color_geno_set_arr.get(i).contains(readinfo_arr.get(index))) { // and
                            color = i;
                            color_geno_set_arr.get(i).add(readinfo_arr.get(index));
                            break;
                        }
                    } else {
                        color_geno_set_arr.add(new HashSet<String>());
                        color_geno_set_arr.get(color_geno_set_arr.size() - 1)
                            .add(readinfo_arr.get(index));

                        color = i;
                        break;
                    }
                }
            }

            if (color > real_max_color) {
                real_max_color = color;
            }

            read_color_arr.set(index, color);
            for (int i = 0; i < readadj_arr.get(index).size(); i++) {
                nb_color_arr.get(readadj_arr.get(index).get(i)).add(color);
            }

            // TODO: [LEFTOVER]
            // System.out.println(max_nb_color);

        }

        // TODO: [LEFTOVER]
        // System.out.println(
        //     "Finished identifying potential full-genome genotypes i.e.: haplotypes.");
        //
        // System.out.println(read_color_arr);

        //
        String null_ref= "*";
        String conf_ref= "";
        for (int i = 1; i < num_loci; i++ ) {
            null_ref = null_ref + "*";
            conf_ref = conf_ref + "?";
        }

        //
        String[] ref_arr = new String[real_max_color+1];
        String[] conf_arr = new String[real_max_color+1];
        for (int i = 0; i <= real_max_color; i++) {
            ref_arr[i] = null_ref;
            conf_arr[i] = conf_ref;
        }

        // TODO: [LEFTOVER]
        // String ss = "1234567";
        // System.out.println(ss.substring(1, ss.length()));

        //
//        for (int i = 0; i < read_pos_2D_arr.size(); i++) {
//        	for (int j = 0; j < read_pos_2D_arr.get(i).size(); j++) {
//        		System.out.print(read_pos_2D_arr.get(i).get(j)+"\t");
//        		
//        	}
//        	System.out.println();
//        }
        
        for (int i = 0; i < read_color_arr.size(); i++) {
            int i_color = read_color_arr.get(i);
            for (int j = 0; j < read_pos_2D_arr.get(i).size(); j++) {
//            	System.out.println(read_pos_2D_arr.get(i).get(j).toString());

                // TODO: [LEFTOVER]
                   // System.out.println(read_pos_2D_arr.get(i).get(j).toString());
            	if (read_pos_2D_arr.get(i).get(j)!=null ) { 
	                int p = read_pos_2D_arr.get(i).get(j);
	                p++;
	                ref_arr[i_color] = ref_arr[i_color].substring(0, (p - 1))
	                    + read_geno_2D_arr.get(i).get(j).toString()
	                    + ref_arr[i_color].substring(p, ref_arr[i_color].length());
            	}
            }

            //
            if (read_pos_2D_arr.get(i).size() > 1) {
                for (int j = 0; j < (read_pos_2D_arr.get(i).size() - 1); j++) {
                    if (Math.abs(read_pos_2D_arr.get(i).get(j + 1) - read_pos_2D_arr.get(i).get(j))
                        == 1) {

                        int p = read_pos_2D_arr.get(i).get(j);
                        p++;
                        try {
                            conf_arr[i_color] = conf_arr[i_color].substring(0, (p - 1))
                                + "-"
                                + conf_arr[i_color].substring(p, conf_arr[i_color].length());

                        } catch (StringIndexOutOfBoundsException e) {
                            System.out.println(conf_arr[i_color].substring(0,(p - 1))
                                + "\t"
                                + p
                                + "\t"
                                + conf_arr[i_color].length());

                        }
                    }
                }
            }

        }

        // TODO: [LEFTOVER]
        // for (int i = 0; i < conf_arr.length; i++) {
        //     System.out.println(ref_arr[i]);
        // }

        this.output_ref_arr = new HashMap<String,Integer>();
        this.conf_ref_arr = new HashMap<String,String>();
        
        int num_window =  this.num_loci/ this.num_loci_window +1;
        if ((this.num_loci % this.num_loci_window )==0) {
        	num_window--;
        }
        ArrayList<ArrayList<String>> ref_reg_2D_arr= new ArrayList<ArrayList<String>>();
        ArrayList<ArrayList<String>> conf_reg_2D_arr= new ArrayList<ArrayList<String>>();
        int completeness_cutoff = this.max_num_gap;
        if (num_window >1 ) {
        	for (int i=0; i <(num_window); i++ ) {
        		ArrayList<String> tmp_ref_arr = new ArrayList<String>();
        		ArrayList<String> tmp_conf_arr = new ArrayList<String>();
        		int index_max = (i+1)* this.num_loci_window;
        		if (((i+1)* this.num_loci_window) > this.num_loci) {
        			index_max = this.num_loci;
        		}
        		
        		for (int j = 0; j <= real_max_color; j++) {
        			if (count(ref_arr[j].substring(this.num_loci_window*i ,index_max )  ) <= completeness_cutoff) {
        				tmp_ref_arr.add( ref_arr[j].substring(this.num_loci_window*i ,index_max )  );
        				tmp_conf_arr.add( conf_arr[j].substring(this.num_loci_window*i,index_max-1 )  );
        			}
        		}
        		ref_reg_2D_arr.add(tmp_ref_arr);
        		conf_reg_2D_arr.add(tmp_conf_arr);
        	}
        	
        	int max_size =0;
        	for (int i=0; i <ref_reg_2D_arr.size(); i++ ) {
        		if (ref_reg_2D_arr.get(i).size()>  max_size) {
        			max_size = ref_reg_2D_arr.get(i).size();
        		}
        	}
        	
        	for (int i=0; i <ref_reg_2D_arr.size(); i++ ) {
        		for (int j=ref_reg_2D_arr.get(i).size();j< max_size;j++) {
        			int max_index = ref_reg_2D_arr.get(i).size();
        			if (max_index>0 ) {
		        	    Random random = new Random(j);
		        	    int s = random.nextInt(max_index)%(max_index+1);
		        	    ref_reg_2D_arr.get(i).add(ref_reg_2D_arr.get(i).get(s));
		        	    conf_reg_2D_arr.get(i).add(conf_reg_2D_arr.get(i).get(s));
        			}else {
        				System.out.println("ERROR: Can not reconstruct haplotype for" + this.vef_file+ "using graph "
        						+ "coloring.\nPlease decrease the Num_Pos_Window or increase Num_Gap_Window!");
        	        	System.exit(0);
        			}
        			
        		}
        	}
        	
//        	for (int i =0; i< ref_reg_2D_arr.size();i++)
//        		System.out.println(ref_reg_2D_arr.get(i));
        	
        	for (int j=0; j< ref_reg_2D_arr.get(0).size(); j++) {
        		String tmp_ref = ref_reg_2D_arr.get(0).get(j); 
        		String tmp_conf =conf_reg_2D_arr.get(0).get(j); 
        		for (int i=1; i <ref_reg_2D_arr.size(); i++ ) {
        			tmp_ref=tmp_ref + ref_reg_2D_arr.get(i).get(j);
        			tmp_conf =tmp_conf+ "?" + conf_reg_2D_arr.get(i).get(j); 
        		}
        		if (this.output_ref_arr.containsKey(tmp_ref)) {
                    this.output_ref_arr.put(tmp_ref, this.output_ref_arr.get(tmp_ref)+1);
                    conf_ref_arr.put(tmp_ref, tmp_conf);

                //
                } else {
                    this.output_ref_arr.put(tmp_ref, 1);
                    this.conf_ref_arr.put(tmp_ref, tmp_conf);
                }
        	}
        	
        	
        }else {
        
	        for (int i = 0; i <= real_max_color; i++) {
	            if (count(ref_arr[i]) <= completeness_cutoff) {	// may implement this in the future
	
	                // TODO: [LEFTOVER]
	                // System.out.println(ref_arr[i]);
	
	                if (this.output_ref_arr.containsKey(ref_arr[i])) {
	                    this.output_ref_arr.put(ref_arr[i], this.output_ref_arr.get(ref_arr[i])+1);
	                    conf_ref_arr.put(ref_arr[i], conf_arr[i]);
	
	                //
	                } else {
	                    this.output_ref_arr.put(ref_arr[i], 1);
	                    this.conf_ref_arr.put(ref_arr[i], conf_arr[i]);
	                }
	
	            } else {
	
	                // TODO: [LEFTOVER]
	                // System.out.println(ref_arr[i] + "\t" + count(ref_arr[i]));
	
	            }
	        }
        }
        
        

        // TODO: [LEFTOVER]
         if (this.conf_ref_arr.isEmpty()) {
        	 String tmp_ref= "1";
        	 String tmp_conf= "";
        	 for (int i=0;i< (this.num_loci-1);i++) {
        		 tmp_ref=tmp_ref+"1";
        		 tmp_conf=tmp_conf+"-";
        	 }
        	 this.conf_ref_arr.put(tmp_ref,  tmp_conf);
        	 this.output_ref_arr.put(tmp_ref, 1);
         }

    }

    public void fileOut(String out_file) throws IOException {
        FileWriter mydata = new FileWriter(out_file,false);
        PrintWriter pw = new PrintWriter(mydata);

        // TODO: [LEFTOVER]
        // System.out.println(output_ref_arr );

        // Changed iterator from Map<K, V> -> K because just getting the key requires less memory.
        for (String entry : this.output_ref_arr.keySet()) {
            String b = this.conf_ref_arr.get(entry);
            String c = "";
            String tmp_str= "";
           
            for (int i = 0; i < b.length(); i++) {
            	tmp_str = entry.substring(i, i + 1);
            	if (tmp_str.equals("*")) {
            		tmp_str="0";
            	}
                c = c + tmp_str + b.substring(i, i + 1); // TODO: (minor) c+=?
            }
            tmp_str =  entry.substring(entry.length() - 1, entry.length()) ;
        	if (tmp_str.equals("*")) {
        		tmp_str="0";
        	}
            c = c+ tmp_str ; // TODO: (minor) c+=?
            pw.write(c + "\t" + this.output_ref_arr.get(entry).toString() + "\n");

            // TODO: [LEFTOVER]
            // System.out.println(c + "\t" + output_ref_arr.get(x).toString());

        }

        pw.flush();
        pw.close();

        // TODO: [LEFTOVER]
        // for (int i = 0; i< real_max_color + 1; i++) {
        //     System.out.println(ref_arr[i]);
        // }

        return;
    }

    /**
     *
     * @return
     */
    public HapConfig hapOut(String[] pool_IDs) {
        //
        int num_global_hap = this.output_ref_arr.size();
        String[][] global_haps_string = new String[num_global_hap][num_loci];
        int[] global_haps_ct = new int[num_global_hap];
        int tot_hap_ct = 0;
        int hap_index = 0;

        //
        for (String entry : this.output_ref_arr.keySet()) {
            String[] var_comp = entry.split("");
            for (int v = 0; v < num_loci; v++) {
                global_haps_string[hap_index][v] = var_comp[v];
            }

            global_haps_ct[hap_index] = this.output_ref_arr.get(entry);
            tot_hap_ct += this.output_ref_arr.get(entry);
            hap_index++;
        }

        //
        double[] global_haps_freq = new double[num_global_hap];
        System.out.println("There are " + num_global_hap + " haplotypes generated by GC.");
        for (int h = 0; h < num_global_hap; h++) {
            global_haps_freq[h] = (double) global_haps_ct[h] / (double) tot_hap_ct;
        }
        return new HapConfig(
            global_haps_string,
            global_haps_freq,
            null,
            this.inpool_site_freqs,
            this.locusInfo,
            this.num_pools,
            null,
            pool_IDs,
            0);

    }


    /**
     *
     * @param gc_hap
     * @return
     */
    int count(String gc_hap) {
        int unknown = 0;
        String[] allele_comp = gc_hap.split("");
        for (String a : allele_comp) {
            if (a.equals("*")) {
                unknown++;
            }
        }

        return unknown;
    }
}
