package PoolHap; 

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ThreadLocalRandom;
import java.util.*;

import PoolHap.Parameters;
import PoolHap.HapConfig;
import PoolHap.SiteInPoolFreqAnno;
import PoolHap.LocusAnnotation;
import PoolHap.RegionEMSolver;

public class DivideConquer {
    /**
     *  @author  Quan Long. Oct 08, 2018
     *
     *  This class divides the full-length haplotype into multiple regions based on linkage support
     *  from reads measured by GC, solves for the regional haplotypes, and merges them. Recursively,
     *  the above procedure will be done multiple times.
     *
     *  There are two strategies on how to divide the genome:
     *  (1) Specify a fixed number of genetic markers (SNPs) for each region
     *  (2) Calculate the number of genetic markers in each region based on estimated LD structure
     *      The LD structure can be estimated by:
     *      (2.1) Another reference panel
     *      (2.2) Information from sequencing reads.
     *  Different layers of divide-and-Conquer may adopt different strategies.
     *
     *  All parameters are in the Object "parameters" to facilitate Machine Learning-based training.
     */

    public Parameters dp;

    public int num_pools;
    public int num_sites;
    public int num_regions_level_I;
    public int[][] regions_level_I;	// num_regsion_level_I x 2 (start and end)
    public int num_regions_level_II;
    public int[][] regions_level_II; // num_regsion_level_II x 2 (start and end)
    
    public int num_regions_level_III;
    public int[][] regions_level_III; // num_regsion_level_III x 2 (start and end)
    
    public int num_regions_level_IV;
    public int[][] regions_level_IV; // num_regsion_level_IV x 2 (start and end)
    
    public int num_regions_level_V ;
    public int[][] regions_level_V; // num_regsion_level_V x 2 (start and end)
    
    public int num_regions_level_VI ;
    public int[][] regions_level_VI; // num_regsion_level_VI x 2 (start and end)
    
    
    public int num_regions_level_VII;
    public int[][] regions_level_VII; // num_regsion_level_VII x 2 (start and end)
    
    public int num_regions_level_VIII;
    public int[][] regions_level_VIII; // num_regsion_level_VIII x 2 (start and end)
    
    public int num_regions_level_final;
    public int[][] regions_level_final; 
    
    public int num_regions_level_final_link;
    public int[][] regions_level_final_link; 
    
    
    
    public int final_level;
    
    
    public String[][] global_haps_gc;
    public double[] global_gc_freq;
    public int num_haps_gc;
    public double[][][] loci_link_freq;
    public HapConfig [] level_I_config;
    public HapConfig [] level_II_config;
    public HapConfig [] level_III_config;
    public HapConfig [] level_IV_config;
    public HapConfig [] level_V_config;
    public HapConfig [] level_VI_config;
    public HapConfig [] level_VII_config;
    public HapConfig [] level_VIII_config;
    
    public HapConfig [] level_final_config;
    public HapConfig [] level_final_link_config;
    
    public int user_set_highest_level;
    
    // Note that we do NOT generate the fields for HapConfigs of each regions. These will be
    // generated on-the-fly during the calculation.

    // The in-pool frequencies and annotations for all sites.
    public SiteInPoolFreqAnno in_pool_sites_freq_anno;
    public String[][] gc_outcome;
    public boolean region_solve; // TODO: remove this when the debugging ends. Quan Long 2019-06-29
    
    public int level_III_IV_region_mismatch_tolerance;
    public int level_V_VI_region_mismatch_tolerance;
    public int level_VII_VIII_region_mismatch_tolerance;
    
    public int[][][] regession_level;  
    public ArrayList<ArrayList<HashMap<Integer,String>>>  regession_level_pos_come_from ;  
    public int[][]  regession_level_num_come_from ;
    double aem_fail_dist_cutoff= 0.005;
    public int aem_fail_max_iter = 5;
    
    

    /**
     *  Constructor that generates regions (i.e.: the dividing plan) based on GC outcomes.
     *
     *
     *  @param frequency_file
     *  @param gc_input_list
     *  @param parameter_file
     *  @param dc_out_file
     */
    
    public DivideConquer(int[][] level_III, int [][] level_IV, 
    		String reg_dc_out_file) {
    	int num_region = level_III.length;
    	int level_1_length= (level_III.length+1)/2;
    	for (int i=0; i<level_III.length; i++ ) {
    		if (Math.pow( 2, i ) >= (double) level_III.length) {
    			num_region= i;
    			break;
    		}
    	}
    	
//    	System.out.println(num_region);
    	this.regession_level = new int[num_region][(level_III.length+1 )/2][2];
    	
    	for (int i=0; i<this.regession_level.length; i++ ) {
    		for (int j=0; j<this.regession_level[i].length; j++ ) {
    			for (int k=0; k<this.regession_level[i][j].length; k++ ) {
    				this.regession_level[i][j][k]=-1;
    			}
    		}
    	}
    	
    	for (int i=0; i< level_III.length /2;i++) {
    		this.regession_level[0][i][0]= level_III[2*i] [0];
    		this.regession_level[0][i][1]= level_III[2*i+1] [1];
    		for (int p=level_III[2*i] [0]; p<=level_III[2*i+1] [1];p++  ) {
    			
    		}
    	}
    	
    	if (level_III.length %2 ==1) {
    		this.regession_level[0][ this.regession_level[0].length-1][0]= 
    				level_III[level_III.length-1] [0];
    		this.regession_level[0][ this.regession_level[0].length-1][1]= 
    				level_III[level_III.length-1] [1];
    	}
    	
    	for (int r=1; r<this.regession_level.length; r++ ) {
    		for (int i=0; i< this.regession_level[r-1].length /2;i++) {
        		this.regession_level[r][i][0]= regession_level[r-1][2*i] [0];
        		this.regession_level[r][i][1]= regession_level[r-1][2*i+1] [1];
        	}
    		int upper_region_len= 0;
        	for (int j=0; j< this.regession_level[r-1].length;j++) {
        		if (this.regession_level[r-1][j][0]!=-1) {
        			upper_region_len++;
        		}
        	}
        	if (upper_region_len %2 ==1) {
        		this.regession_level[r][upper_region_len/2][0]= 
        				this.regession_level[r-1][upper_region_len-1] [0];
        		this.regession_level[r][upper_region_len/2][1]= 
        				this.regession_level[r-1][upper_region_len-1] [1];
        	}
    	}
    	
    	for (int i=0; i<this.regession_level.length; i++ ) {
    		for (int j=0; j<this.regession_level[i].length; j++ ) {
    			System.out.print(this.regession_level[i][j][0]+";"+ this.regession_level[i][j][1]+"\t");
    		}
    		System.out.println();
    	}
    }
    
    public DivideConquer(int[][] level_III, int [][] level_IV, 
    		String reg_dc_out_file, int num_regions_merge, 
    		HashMap<Integer, Integer> index_var_pos_dict) {
    	int num_level = level_III.length;
    	for (int i=0; i<level_III.length; i++ ) {
    		if (Math.pow( num_regions_merge , i ) >= level_III.length) {
    			num_level= i;
    			break;
    		}
    	}
    	if( num_level==0) {
    		num_level++;
    	}
//    	System.out.println(num_region);
    	if ((level_III.length %num_regions_merge) ==0 ) {
    		this.regession_level = new int[num_level]
    				[(level_III.length+1 )/num_regions_merge ][2];
    		this.regession_level_num_come_from= new int[num_level]
    				[(level_III.length+1 )/num_regions_merge ];
    	} else {
    		this.regession_level = new int[num_level]
    				[(level_III.length)/num_regions_merge +1][2];
    		this.regession_level_num_come_from = new int[num_level]
    				[(level_III.length)/num_regions_merge +1];
    	}
    	
    	this.regession_level_pos_come_from= new ArrayList<ArrayList<HashMap<Integer,String>>>();
    	for (int i=0;i<this.regession_level.length; i++) {
    		this.regession_level_pos_come_from .add(new ArrayList<HashMap<Integer,String>>());
    		for (int j=0;j<this.regession_level[i].length; j++) {
    			this.regession_level_pos_come_from.get(i).add(new HashMap<Integer,String> ());
    		}
    	}
    	
    	for (int i=0; i<this.regession_level.length; i++ ) {
    		for (int j=0; j<this.regession_level[i].length; j++ ) {
    			this.regession_level_num_come_from[i][j]=-1;
    			for (int k=0; k<this.regession_level[i][j].length; k++ ) {
    				this.regession_level[i][j][k]=-1;
    				
    			}
    		}
    	}
    	
    	for (int i=0; i< level_III.length /num_regions_merge;i++) {
    		this.regession_level[0][i][0]= level_III[num_regions_merge*i] [0];
    		this.regession_level[0][i][1]= level_III
    				[num_regions_merge*i+num_regions_merge-1] [1];
    		for (int j=0; j< (num_regions_merge); j++) {
    			for (int p =level_III[i*num_regions_merge+ j ][0]; 
    					p<=  level_III[i*num_regions_merge+ j ][1]; p++) {
//    				if (index_var_pos_dict.get(p)==null) {
//        				System.out.print("***1***  "+p);
//        			}
    				this.regession_level_pos_come_from.get(0).get(i).put(index_var_pos_dict.get(p) , 
    						Integer.toString(i*num_regions_merge+ j)+ ":"+ Integer.toString(j));
    			}
    		}
    	}
    	
    	
    	if (level_III.length %num_regions_merge !=0 ) {
    		this.regession_level[0][ this.regession_level[0].length-1][0]= 
    				level_III[level_III.length-level_III.length% num_regions_merge] [0]; 
    		this.regession_level[0][ this.regession_level[0].length-1][1]= 
    				level_III[level_III.length-1] [1];
//    		for (int j= level_III.length-1 ; j>level_III.length-1-
//    				level_III.length %num_regions_merge  ; j-- ) {
    		for (int j=level_III.length-level_III.length% num_regions_merge ;
    				j<= (level_III.length-1) ; j++ ) {
    			for (int p =level_III[j ][0]; 
    					p<=  level_III[j ][1]; p++) {
    				if (index_var_pos_dict.get(p)==null) {
//        				System.out.print("******2********  "+p);
        			}
    				this.regession_level_pos_come_from.get(0).
    				get(this.regession_level[0].length-1).put(index_var_pos_dict.get(p) , 
    						Integer.toString(j)+ ":"+Integer.toString(j-level_III.length+
    	    				level_III.length %num_regions_merge)  );
    			}
    		}
    	}
    	
    	for (int r=1; r<this.regession_level.length; r++ ) {
    		int upper_region_len= 0;
        	for (int j=0; j< this.regession_level[r-1].length;j++) {
        		if (this.regession_level[r-1][j][0]!=-1) {
        			upper_region_len++;
        		}
        	}
    		for (int i=0; i< upper_region_len /num_regions_merge;i++) {
        		this.regession_level[r][i][0]= this.regession_level[r-1]
        				[num_regions_merge*i] [0];
        		this.regession_level[r][i][1]= this.regession_level[r-1]
        				[num_regions_merge*i+num_regions_merge-1] [1];
        		for (int j= num_regions_merge*i; j<= (num_regions_merge*i+num_regions_merge-1);j++) {
        			for (int p =this.regession_level[r-1][j ][0]; 
        					p<=  this.regession_level[r-1][j ][1]; p++) {
//        				if (index_var_pos_dict.get(p)==null) {
//	        				System.out.println("*****3*********  "+p+"&&"+r+"$$"+i+"??");
//	        			}
        				this.regession_level_pos_come_from.get(r).
        				get(i).put(index_var_pos_dict.get(p) ,
        						Integer.toString(j)+":"+Integer.toString(j - num_regions_merge*i));
        			}
        		}
        	}
    		
        	if ( (upper_region_len % num_regions_merge) !=0) {
	        	this.regession_level[r][upper_region_len/num_regions_merge][0]= 
	        			this.regession_level[r-1][upper_region_len-
	        				upper_region_len % num_regions_merge ] [0];
	        	this.regession_level[r][upper_region_len/num_regions_merge][1]= 
	        				this.regession_level[r-1][upper_region_len-1] [1];
	        	for (int j=upper_region_len-upper_region_len % num_regions_merge; 
	        			j<= (upper_region_len-1);j++) {
	        		for (int p =this.regession_level[r-1][j ][0]; 
        					p<= this.regession_level[r-1][j ][1]; p++) {
//	        			if (index_var_pos_dict.get(p)==null) {
//	        				System.out.print("*******4*******  "+p);
//	        			}
	        			this.regession_level_pos_come_from.get(r).
        				get(upper_region_len/num_regions_merge).
        				put(index_var_pos_dict.get(p) , 
        					Integer.toString(j)+":"+
        					Integer.toString(j-upper_region_len+upper_region_len % num_regions_merge ));
	        		}
	        	}
        	}
    	}
    	
    	for (int i=0; i<this.regession_level.length; i++ ) {
    		for (int j=0; j<this.regession_level[i].length; j++ ) {
    			System.out.print(this.regession_level[i][j][0]+";"+ this.regession_level[i][j][1]+"\t");
    		}
    		System.out.println();
    		
    	}
    	
    	for (int i=0; i<this.regession_level.length; i++ ) {
    		for (int j=0; j<this.regession_level[i].length; j++ ) {
    			HashSet<String>  region_set = new  HashSet<String>();
    			for (Integer pos : this.regession_level_pos_come_from.get(i).get(j).keySet()) {
//    				System.out.print(pos+";"+  this.regession_level_pos_come_from.get(i).get(j).get(pos)+"\t");
    				region_set.add(this.regession_level_pos_come_from.get(i).get(j).get(pos)  );
    			}
    			this.regession_level_num_come_from[i][j]= region_set.size();
    		}
//    		System.out.println();
    	}
    	
    	for (int i=0; i<this.regession_level_num_come_from.length; i++ ) {
    		for (int j=0; j<this.regession_level_num_come_from[i].length; j++ ) {
    			System.out.print(this.regession_level_num_come_from[i][j]+"\t");
    		}
    		System.out.println();
    		
    	}
    	
    	
    	for (int i=0; i<this.regession_level.length; i++ ) {
    		for (int j=0; j<this.regession_level[i].length; j++ ) {
    			for (Integer pos : this.regession_level_pos_come_from.get(i).get(j).keySet()) {
    				System.out.print(pos+";"+  this.regession_level_pos_come_from.get(i).get(j).get(pos)+"\t");
    			}
    		}
    		System.out.println();
    	}
    	
    }
	
    public DivideConquer(
        String frequency_file,
        String[] gc_input_list,
        String parameter_file,
        String dc_out_file) {
    	

        try {
            /*
             *  Load files.
             */
        	        	
        	int num_freq_sites= 0;
        	
        	BufferedReader br = new BufferedReader(new FileReader(frequency_file));
            String line = br.readLine();

            while (line.startsWith("#")) {
                line = br.readLine();
            }
            
            line = br.readLine();

            while (line != null) {
            	num_freq_sites++;
                line = br.readLine();
            }

            br.close();
            
            this.dp = new Parameters(parameter_file); 
            
            if (dp.max_level_I_region_size >= num_freq_sites) {
            	dp.max_level_I_region_size = num_freq_sites-1;
            	dp.min_level_I_region_size = num_freq_sites-1;
            	dp.max_level_II_region_size = num_freq_sites-1;
            	dp.min_level_II_region_size = num_freq_sites-1;
            }
            
            System.out.println("Finished loading the PoolHapX parameter file from " + parameter_file);
            //load_gc_outcome(parse_gc_input(gc_input_list)); 
            // above was removed and replaced by the line below by Quan Long 2019-07-07
            this.load_gc_outcome((gc_input_list));
            System.out.println("Finished loading graph-coloring files from " + dp.inter_dir+"gcf");
            System.out.println("Number of pools = "
                + this.num_pools
                + "\tNumber of segregating sites = "
                + this.num_sites);

            /*
             *  Division?
             */
            
           
            
            this.generate_dividing_plan_two_level(this.gc_outcome);
            if (this.dp.aem_maximum_level>=3) {
	            this.generate_dividing_plan_level_III();
	            this.generate_dividing_plan_level_IV();
            }
            if (this.dp.aem_maximum_level>=5) {
	            this.generate_dividing_plan_level_V();
	            this.generate_dividing_plan_level_VI();
            }
            if (this.dp.aem_maximum_level>=7) {
	            this.generate_dividing_plan_level_VII();
	            this.generate_dividing_plan_level_VIII();
            }
            
            System.out.println(  "----------------------");
            System.out.println(  this.num_regions_level_I);
            System.out.println(  this.num_regions_level_III);
            System.out.println(  this.num_regions_level_V);
            System.out.println(  this.num_regions_level_VII);
            System.out.println(  "-----------------------");
        	this.final_level =7;
            if (this.num_regions_level_III ==0) {
            	this.final_level =1;
            } else if (this.num_regions_level_V ==0) {
            	this.final_level =3;
            } else if (this.num_regions_level_VII ==0) {
            	this.final_level =5;
            }
            
            if (this.final_level > this.dp.aem_maximum_level) {
            	this.final_level = this.dp.aem_maximum_level;
            }
            
            if (this.final_level ==7) {
            	this.num_regions_level_final= this.num_regions_level_VII;
            	this.regions_level_final= this.regions_level_VII.clone();
            	this.num_regions_level_final_link= this.num_regions_level_VIII;
            	this.regions_level_final_link= this.regions_level_VIII.clone();
            }else if (this.final_level ==5) {
            	this.num_regions_level_final= this.num_regions_level_V;
            	this.regions_level_final= this.regions_level_V.clone();
            	this.num_regions_level_final_link= this.num_regions_level_VI;
            	this.regions_level_final_link= this.regions_level_VI.clone();
            } else if (this.final_level ==3) {
            	this.num_regions_level_final= this.num_regions_level_III;
            	this.regions_level_final= this.regions_level_III.clone();
            	this.num_regions_level_final_link= this.num_regions_level_IV;
            	this.regions_level_final_link= this.regions_level_IV.clone();
            } else if (this.final_level ==1) {
            	this.num_regions_level_final= this.num_regions_level_I;
            	this.regions_level_final= this.regions_level_I.clone();
            	this.num_regions_level_final_link= this.num_regions_level_II;
            	this.regions_level_final_link= this.regions_level_II.clone();
            }
            
          
            this.loci_link_freq=new double [num_pools][num_sites-1][4];
            // TODO: LEFTOVER ML 20190702
            // // If can be divided...
            // if (this.region_solve) {

            this.output_current_DC_plan(dc_out_file);
            System.out.println("The current dividing plan: ");
            System.out.println(
                "\nFinished generating divide-and-conquer plan. The plan has been written to "
                + dc_out_file);

            // TODO: LEFTOVER ML 20190702
            // // If not enough gaps...
            // }
            // else {
            //     System.out.println(
            //         "There are not enough locations of linkage uncertainty to run regional "
            //         + "reconstruction. Skipping straight to linkage-informed LASSO regression.");
            // }

            this.in_pool_sites_freq_anno = new SiteInPoolFreqAnno(frequency_file);

        } catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    public DivideConquer() {
		// TODO Auto-generated constructor stub
	}

	/**
     * @author Chen Cao 2019-09
     * Level I and Level II have provided regional haplotypes for AEM.
     * link_regional_AEM() function link the level I regional haplotypes
     * according the level II.
     * After that, AEM is used for the linked larger level I regions.
     * 
     */  
    public void generate_dividing_plan_level_III() {

	    	this.num_regions_level_III= (this.num_regions_level_I +1)/2; 
	    	if (this.num_regions_level_I ==1) {
	    		this.num_regions_level_III =0;
	    	}
	    	if (this.num_regions_level_III!=0 ) {
		    	this.regions_level_III = new int[num_regions_level_III][2];
		    	for (int i=0; i< this.num_regions_level_I /2;i++) {
		    		this.regions_level_III[i][0]= regions_level_I[2*i] [0];
		    		this.regions_level_III[i][1]= regions_level_I[2*i+1] [1];
		    	}
		    	if (this.num_regions_level_I %2 ==1) {
		    		this.regions_level_III[ this.num_regions_level_III-1][0]= 
		    				regions_level_I[this.num_regions_level_I-1] [0];
		    		this.regions_level_III[ this.num_regions_level_III-1][1]= 
		    				regions_level_I[this.num_regions_level_I-1] [1];
		    	}
	    	}
    }
    
    public void generate_dividing_plan_level_IV() { 
		    	this.num_regions_level_IV= (this.num_regions_level_I )/2;
		    	this.regions_level_IV = new int[num_regions_level_IV][2];
		    	for (int i=0; i< (this.num_regions_level_I-1 )/2;i++) {
		    		this.regions_level_IV[i][0]= regions_level_I[2*i+1] [0];
		    		this.regions_level_IV[i][1]= regions_level_I[2*i+2] [1];
		    	}
		    	if ((this.num_regions_level_I) %2 ==0) {
		    		this.regions_level_IV[ this.num_regions_level_IV-1][0]= 
		    				regions_level_I[this.num_regions_level_I-1] [0];
		    		this.regions_level_IV[ this.num_regions_level_IV-1][1]= 
		    				regions_level_I[this.num_regions_level_I-1] [1];
		    	}
    }
    
    public void generate_dividing_plan_level_V() {
    	
	    	this.num_regions_level_V= (this.num_regions_level_III +1)/2;
	    	if (this.num_regions_level_III ==1) {
	    		this.num_regions_level_V =0;
	    	}
	    	if (this.num_regions_level_V!=0 ) {
		    	this.regions_level_V = new int[num_regions_level_V][2];
		    	for (int i=0; i< this.num_regions_level_III /2;i++) {
		    		this.regions_level_V[i][0]= regions_level_III[2*i] [0];
		    		this.regions_level_V[i][1]= regions_level_III[2*i+1] [1];
		    	}
		    	if (this.num_regions_level_III %2 ==1) {
		    		this.regions_level_V[ this.num_regions_level_V-1][0]= 
		    				regions_level_III[this.num_regions_level_III-1] [0];
		    		this.regions_level_V[ this.num_regions_level_V-1][1]= 
		    				regions_level_III[this.num_regions_level_III-1] [1];
		    	}
	    	}
    }
    
    public void generate_dividing_plan_level_VI() {
    	if (this.num_regions_level_III> 1){
	    	this.num_regions_level_VI= (this.num_regions_level_III )/2;
	    	this.regions_level_VI = new int[num_regions_level_VI][2];
	    	for (int i=0; i< (this.num_regions_level_III-1 )/2;i++) {
	    		this.regions_level_VI[i][0]= regions_level_III[2*i+1] [0];
	    		this.regions_level_VI[i][1]= regions_level_III[2*i+2] [1];
	    	}
	    	
	    	if ((this.num_regions_level_III) %2 ==0) {
	    		this.regions_level_VI[ this.num_regions_level_VI-1][0]= 
	    				regions_level_III[this.num_regions_level_III-1] [0];
	    		this.regions_level_VI[ this.num_regions_level_VI-1][1]= 
	    				regions_level_III[this.num_regions_level_III-1] [1];
	    	}
    	}
    }
    
    
    
    public void generate_dividing_plan_level_VII() {
    	
    	this.num_regions_level_VII= (this.num_regions_level_V +1)/2;
    	if (this.num_regions_level_V ==1) {
    		this.num_regions_level_VII =0;
    	}
    	if (this.num_regions_level_VII!=0 ) {
	    	this.regions_level_VII = new int[num_regions_level_VII][2];
	    	for (int i=0; i< this.num_regions_level_V /2;i++) {
	    		this.regions_level_VII[i][0]= regions_level_V[2*i] [0];
	    		this.regions_level_VII[i][1]= regions_level_V[2*i+1] [1];
	    	}
	    	if (this.num_regions_level_V %2 ==1) {
	    		this.regions_level_VII[ this.num_regions_level_VII-1][0]= 
	    				regions_level_V[this.num_regions_level_V-1] [0];
	    		this.regions_level_VII[ this.num_regions_level_VII-1][1]= 
	    				regions_level_V[this.num_regions_level_V-1] [1];
	    	}
    	}
    	
    }

	public void generate_dividing_plan_level_VIII() {
		if (this.num_regions_level_V> 1){
	    	this.num_regions_level_VIII= (this.num_regions_level_V )/2;
	    	this.regions_level_VIII = new int[num_regions_level_VIII][2];
	    	for (int i=0; i< (this.num_regions_level_V-1 )/2;i++) {
	    		this.regions_level_VIII[i][0]= regions_level_V[2*i+1] [0];
	    		this.regions_level_VIII[i][1]= regions_level_V[2*i+2] [1];
	    	}
	    	
	    	if ((this.num_regions_level_V) %2 ==0) {
	    		this.regions_level_VIII[ this.num_regions_level_VIII-1][0]= 
	    				regions_level_V[this.num_regions_level_V-1] [0];
	    		this.regions_level_VIII[ this.num_regions_level_VIII-1][1]= 
	    				regions_level_V[this.num_regions_level_V-1] [1];
	    	}
		}
	}
    
    public void link_regional_AEM_level_I(){
    
   
    }

    /**
     *  Read files recording the outcome of GC:
     *
     *  0?0-0-0?1-0-0-0-0-0-0-1-0-0-0-0-0?0-0-0?1?0?0-0-0-1-0-0?0-0\t25
     *  0?0-0-0?1-0-0-0-0-0-0-1-0-1-1-0-1?0-0-0?1?0?0-0-0-0-0-0?1-1\t10
     *  1?1-0-1?1-0-0-0-0-0-0-1-0-1-1-0-1?0-1-0?0?0?0-1-0-0-0-0?1-1\t8
     *
     *  It will return GC outcomes in the form of strings.
     *
     *  @param graph_coloring_outcome_files
     */
    
    
    
    
    public void load_gc_outcome(String[] graph_coloring_outcome_files) {
        this.num_pools = graph_coloring_outcome_files.length;
        this.gc_outcome = new String[num_pools][];
        try {
            for (int p_index = 0; p_index < num_pools; p_index++) {
                ArrayList<String> haps_list = new ArrayList<String>();
                BufferedReader br = new BufferedReader(
                    new FileReader(new File(graph_coloring_outcome_files[p_index])));

                String line = br.readLine();

                // TODO: [ReconEP]:: extract to static utils method.
                while (line.startsWith("#")) {
                    line = br.readLine();
                }

                // Changed from original code (+1 -> -1) because at the end of the draft haplotype
                // variant.
                this.num_sites = (line.split("\t")[0].length() + 1) / 2;

                // Composition, the output_ref_arr.get(x), which seems to be the count of that
                // draft.
                while (line != null) {
                    // Haplotype in the pool, is also outputted. Need to NOT read this into the
                    // variant composition.
                    haps_list.add(line);
                    line = br.readLine();
                }

                this.gc_outcome[p_index] = new String[haps_list.size()];
                for (int h = 0; h < haps_list.size(); h++) {
                    this.gc_outcome[p_index][h] = haps_list.get(h);
                }

                br.close();
            }

        } catch(Exception e) {
            e.printStackTrace();
        }
    }


    /**
     *
     *  @param gc_pathes_file
     *  @return
     */
    public String[] parse_gc_input(String gc_pathes_file) {
        ArrayList<String> pathes = new ArrayList<String>();
        try{
            BufferedReader br = new BufferedReader(new FileReader(gc_pathes_file));
            String line = br.readLine();

            // TODO: [ReconEP]:: extract to static utils method.
            while (line.startsWith("#")) { // skip headers, if any
                line = br.readLine();
            }

            while (line != null) {
                pathes.add(line);
                line = br.readLine();
            }

            this.num_pools = pathes.size();
            br.close();

        } catch(Exception e) {
            e.printStackTrace();
        }

        String[] gc_input_files = new String[this.num_pools];
        for (int p = 0; p < this.num_pools; p++) {
            gc_input_files[p] = pathes.get(p);
        }
        return gc_input_files;
    }


    /**
     *
     * @param graph_coloring_outcome
     */
    public void generate_dividing_plan_two_level(String[][] graph_coloring_outcome){
        // Step 1) Generate list of positions where there are linkage uncertainties i.e.: gaps.
        int[] gap_positions = identify_gaps(graph_coloring_outcome);

        // Step 2) Make the level 1 regions i.e.: windows of variants based on the generated gaps.
        // These are the boundaries of the regions/windows of variants.
        // TODO: [ReconEP]:: again, each step should be helpers?
        ArrayList<Integer> region_cuts = new ArrayList<Integer>();

        // This is the length of (number of variants in) the region/windows of variants so far.
        
        
        int curr_unmatched_len=0;
        
        //Make sure that PoolHapX at least has two level I regions.
//        if (dp.max_level_I_region_size >=  this.num_sites) {
//        	dp.max_level_I_region_size= (this.num_sites+1)/2; 
//        	dp.min_level_I_region_size= (this.num_sites+1)/2; 
//        }
//        
        for (int gap = 0; gap < gap_positions.length; gap++) { // For each gap...
            // The start position of the region currently being made.
            int previous_gap_position = (gap == 0)?0:gap_positions[gap - 1]; // No linting?

            // The number of variants in the region currently being made.
            int gap_len = gap_positions[gap] - previous_gap_position;
            // regardless whether the above is executed, add the gap_len to curr_unmatched_len
            curr_unmatched_len += gap_len;
            // If the newly updated curr_unmatched_len contains enough variants to fall within the
            // acceptable range, form the next region.
            if (inbetween(curr_unmatched_len, 
                dp.min_level_I_region_size,  dp.max_level_I_region_size)) {
                // Create the region.
                curr_unmatched_len = 0;
                region_cuts.add(gap_positions[gap]);
            } // ...otherwise, if the region currently being made is larger than the max allowable
            // size, split it up so that everything after the max size is added to the next region.
            else if (curr_unmatched_len > dp.max_level_I_region_size) {
                curr_unmatched_len = curr_unmatched_len - dp.max_level_I_region_size;
                int last_cut= (region_cuts.size()==0)?0:(region_cuts.get(region_cuts.size()-1));
                region_cuts.add(last_cut + dp.max_level_I_region_size);
                // the current_unmatched_len may still be larger than max_level_I_region_size
                // because of the gap_len just added is very large. 
                while(curr_unmatched_len > dp.max_level_I_region_size) {
                    curr_unmatched_len = curr_unmatched_len - dp.max_level_I_region_size;
                    last_cut= region_cuts.get(region_cuts.size()-1);
                    region_cuts.add(last_cut + dp.max_level_I_region_size);
                }// If the remained curr_unmatched_len contains enough variants to fall within the
                // acceptable range, form the next region.
                if (inbetween(curr_unmatched_len, 
                    dp.min_level_I_region_size,  dp.max_level_I_region_size)) {
                    // Create the region.
                    last_cut= region_cuts.get(region_cuts.size()-1);
                    region_cuts.add(last_cut+curr_unmatched_len);
                    curr_unmatched_len = 0;
                }
            } else { // curr_unmatched_len < parameters.min_level_I_region_size
                // do NOT form a region;
            }
        }

        this.num_regions_level_I = region_cuts.size() + 1;

        // The start and end positions of each region.
        this.regions_level_I = new int[num_regions_level_I][2];
        this.regions_level_I[0][0] = 0;

        // Note that because there are max region sizes, not every region will start and/or end at a
        // gap position.

        // Gap positions are merely guidelines to determine linkage regions.
        for (int r = 0; r < num_regions_level_I - 1; r++) {
            this.regions_level_I[r][1] = region_cuts.get(r) - 1;
            this.regions_level_I[r + 1][0] = region_cuts.get(r);
        }

        this.regions_level_I[num_regions_level_I - 1][1] = this.num_sites - 1;

        // Step 3) Make the level 2 regions i.e.: windows of variants that overlap with the end of
        // the ith level 1 window and the start of the (i + 1)th level 1 window.
        this.num_regions_level_II = this.num_regions_level_I - 1;
        if (this.num_regions_level_II!=0) {
	        this.regions_level_II = new int[num_regions_level_II][2];
	
	        // !!! This is the midpoint of the 0th level 1 region.
	        this.regions_level_II[0][0] = (this.regions_level_I[0][1] - this.regions_level_I[0][0]) / 2;
	        for (int r = 0; r < num_regions_level_II - 1; r++) {
	            int tmp_end = (this.regions_level_I[r+1][0] + this.regions_level_I[r + 1][1]) / 2;
	            int tmp_reg_2_len = tmp_end - this.regions_level_II[r][0] + 1;
	            if (inbetween(
	                tmp_reg_2_len,
	                dp.min_level_II_region_size,
	                dp.max_level_II_region_size)) {
	
	                this.regions_level_II[r][1] = tmp_end;
	
	            } else if (tmp_reg_2_len > dp.max_level_II_region_size) {
	                this.regions_level_II[r][1] = this.regions_level_II[r][0]
	                    + dp.max_level_II_region_size;
	
	            } else {
	                this.regions_level_II[r][1] = tmp_end;
	                while (tmp_reg_2_len < dp.min_level_II_region_size) {
	                    this.regions_level_II[r][1]++;
	                    tmp_reg_2_len = this.regions_level_II[r][1] - this.regions_level_II[r][0] + 1;
	
	                }
	            }
	
	            this.regions_level_II[r + 1][0] = this.regions_level_II[r][1] + 1;
	        }
	
	        this.regions_level_II[num_regions_level_II - 1][1] = 
	        		(this.regions_level_I[num_regions_level_II][0] + 
	        				this.regions_level_I[num_regions_level_II][1]) / 2;
        }
    }


    /**
     *
     *  @param value
     *  @param min
     *  @param max
     *  @return
     */
    public boolean inbetween(int value, int min, int max) {
        return (value >= min && value <= max);
    }
    

    public HashMap<Integer, String>   gs_map(String gs_var_pos) throws IOException {
    	BufferedReader br = new BufferedReader(new FileReader(gs_var_pos));
        String currLine = br.readLine(); // skip header
        currLine = br.readLine();
        int pos_index=0 ; 
        // Read each line into a dictionary.
        HashMap<Integer, String>pos_dict  =new  HashMap<Integer, String>();
        while (currLine != null) {
            pos_dict.put(pos_index, currLine.split("\t")[0]);

            // TODO: [LEFTOVER]
            // System.out.println(currLine.split(";")[1] + "\t" +  pos_index);

            pos_index++; // move to next variant position index
            currLine = br.readLine(); // read next line
        }
        br.close();
    	return pos_dict;
    	
    }
    
    public HashMap<Integer, Integer>   gs_map_pos(String gs_var_pos) throws IOException {
    	BufferedReader br = new BufferedReader(new FileReader(gs_var_pos));
        String currLine = br.readLine(); // skip header
        currLine = br.readLine();
        int pos_index=0 ; 
        // Read each line into a dictionary.
        HashMap<Integer, Integer>pos_dict  =new  HashMap<Integer, Integer>();
        while (currLine != null) {
            pos_dict.put(pos_index,
            		Integer.parseInt(currLine.split(";")[1]));

            // TODO: [LEFTOVER]
            // System.out.println(currLine.split(";")[1] + "\t" +  pos_index);

            pos_index++; // move to next variant position index
            currLine = br.readLine(); // read next line
        }
        br.close();
    	return pos_dict;
    	
    }

    /** TODO: remove this funtion when everything is stablized. 
     *  Identify gaps from the GC outcome.
     *  The GC outcome file is composed of inferred haplotypes in the format of:
     *  		0-1-0-1-0?1-0-1\t5
     *  where 0/1 stands for alleles, "-" stands for linked, and "?" stands for gap. The number
     *  after \t stands for the count
     *
     *  There are two criteria:
     *      (1) the gap in a single pool, largely due to the sequencing gap
     *      (2) the shared gap in all the pools, largely due to the long physical distance between
     *          two sites.
     *
     *  Parameters involved:
     *      double gap_all_pool_cutoff
     *      double gap_inpool_cutoff
     *
     *  @param graph_coloring_outcome
     *  @return
     */
    /**
    public int[] identify_gaps_Lauren_iterative_tobe_deleted_later(String[][] graph_coloring_outcome) {
        ArrayList<Integer> gap_indexes = new ArrayList<>();

        // Step 1) Count up the number of gaps within each and between all of the raw GC haplotypes
        // and all of the pools.
        // TODO: [ReconEP]:: separate helpers for each step?
        // The count at each potential gap in each pool. "this.num_sites-1" potential gaps.
        int[][] gap_counts_inpool = new int[this.num_pools][this.num_sites - 1];

        // The cumulative count at each potential gap position.
        int[] gap_counts_all = new int[this.num_sites - 1];

        // The number of types of raw GC haplotypes in each pool.
        double[] num_haps_inpool = new double[this.num_pools];
        double num_haps_all = 0.0;
        HashMap<String,Integer> hap_tracker = new HashMap<String,Integer>();
        int tot_ct = 0;
        for (int p = 0; p < this.num_pools; p++) { // for each pool...
            String[] haps = graph_coloring_outcome[p]; // the list of raw GC haplotypes
            num_haps_all += haps.length;
            num_haps_inpool[p] = haps.length;
            for (int h = 0; h < haps.length; h++){ // ...for each raw GC haplotype...
                String curr_vc = "";
                int hap_ct = Integer.parseInt(haps[h].split("\t")[1]);
                for(int k = 0; k < this.num_sites - 1; k++) { // ...for each potential gap...
                    curr_vc += haps[h].charAt(k * 2);

                    // If the linkage between the two variant positions is uncertain
                    // i.e.: gap is present...
                    if (haps[h].charAt(k * 2 + 1) == '?'){
                        gap_counts_inpool[p][k] += hap_ct; // increment the count within-pool and globally
                        gap_counts_all[k] += hap_ct;
                    }
                }

                curr_vc += haps[h].charAt((this.num_sites - 1) * 2);
                if (!hap_tracker.containsKey(curr_vc)) {
                    hap_tracker.put(curr_vc, hap_ct);

                } else {
                    int new_ct = hap_tracker.get(curr_vc) + hap_ct;
                    hap_tracker.put(curr_vc, new_ct);
                }

                tot_ct += hap_ct;
            }
        }

        // Step 2) Check if the potential gaps meet the within- OR between-pool frequency
        // thresholds.
        HashSet<Integer> gap_indexes_set = new HashSet<Integer>();
        double avg_region_size = this.num_sites;
        double across_cutoff = 1;
        double local_cutoff = 1;

        // While there aren't enough gaps to divide up the GC haplotypes...
        while (across_cutoff >= this.dp.gap_all_pool_cutoff
            || local_cutoff >= this.dp.gap_inpool_cutoff) {

            for (int k = 0; k <this.num_sites - 1; k++) {
                if((double) gap_counts_all[k] / num_haps_all >= across_cutoff
                    && !gap_indexes_set.contains(k)) {

                    gap_indexes.add(k);
                    gap_indexes_set.add(k);
                }

                for (int p = 0; p < this.num_pools; p++) {
                    if ((double) gap_counts_inpool[p][k] / num_haps_inpool[p] >= local_cutoff
                        && !gap_indexes_set.contains(k)){

                        gap_indexes.add(k);
                        gap_indexes_set.add(k);
                    }
                }
            }

            avg_region_size = (double) this.num_sites / (double) (gap_indexes.size() + 1);
            if (avg_region_size > this.dp.max_level_I_region_size) {
                break;
            }

            // Relax the cutoffs so that there many be more gaps.
            across_cutoff -= this.dp.gap_support_step;

            // The local cutoff is identical in case there are rare haplotypes with uncommon linkage
            // patterns.
            local_cutoff -= this.dp.gap_support_step;
        }

        System.out.println("The final across-pool cutoff is "
            + across_cutoff
            + " and the within-pool cutoff is "
            + local_cutoff
            + ". The average region size is "
            + avg_region_size
            + " segregating sites long.");

        int[] gap_indexes_array;

        // If there are enough gaps to run regional reconstruction....
        if (avg_region_size <= this.dp.max_level_I_region_size) {
            gap_indexes_array = new int[gap_indexes.size()];

            // Add the last site index as the final "gap" so that all the sites are within gaps.
            gap_indexes.add(this.num_sites - 1);

            // This is not useful for the moment, but keep the data integrity for potential future
            // use.
            gap_indexes_set.add(this.num_sites - 1);

            // Clean the outcome and return.
            for (int i = 0; i < gap_indexes.size(); i++) {
                gap_indexes_array[i] = gap_indexes.get(i);
            }
            Arrays.sort(gap_indexes_array); // sort the indices

            this.global_haps_gc = new String[hap_tracker.size()][this.num_sites];
            this.global_gc_freq = new double[hap_tracker.size()];
            this.num_haps_gc = hap_tracker.size();
            int hap_index = 0;
            for (String curr_vc : hap_tracker.keySet()) {
                String[] tmp = curr_vc.split("");
                for (int l = 0; l < this.num_sites; l++) {
                    this.global_haps_gc[hap_index][l] = tmp[l];
                }

                this.global_gc_freq[hap_index] = (double) hap_tracker.get(curr_vc)
                    / (double) tot_ct;

                hap_index++;
            }

            this.region_solve = true;

        } else { // ...otherwise, we can just reduce the candidates using LASSO
            gap_indexes_array = new int[1];
            gap_indexes_array[0] = -1;
            this.region_solve = false;
        }

        return gap_indexes_array;
    }
    
    
    
    */

    /*
     *  Identify gaps from the GC outcome.
     *  The GC outcome file is composed of inferred haplotypes in the format of:
     * 		0-1-0-1-0?1-0-1\t30
     *  where 0/1 stands for alleles, "-" stands for linked, and "?" stands for gap. The number
     *  after \t stands for the count of the hap
     *
     *  There are two criteria:
     * 	 (1) the gap in a single pool, largely due to the sequencing gap
     * 	 (2) the shared gap in all the pools, largely due to the long physical distance between two
     *       sites.
     *
     *  Parameters involved:
     * 	 double gap_all_pool_cutoff
     * 	 double gap_inpool_cutoff
     */
    
    
    public int[] identify_gaps(String[][] graph_coloring_outcome) {
        ArrayList<Integer> gap_indexes = new ArrayList<>();
        // Step 1) Count up the number of gaps within each and between all of the raw GC haplotypes
        // and all of the pools.
        // The count at each potential gap in each pool. "this.num_sites-1" potential gaps.
        
        int[][] gap_counts_inpool = new int[this.num_pools][this.num_sites-1];
        // The cumulative count at each potential gap position.
        int[] gap_counts_all=new int[this.num_sites-1];
        // The number of types of raw GC haplotypes in each pool.
        double[] num_haps_inpool=new double[this.num_pools];
  
        double num_haps_all = 0.0;
        HashMap<String,Integer> hap_tracker = new HashMap<String,Integer>();
        int tot_ct = 0;
        for (int p = 0; p < this.num_pools; p++) { // for each pool...
            String[] haps = graph_coloring_outcome[p]; // the list of raw GC haplotypes
            num_haps_all += haps.length;
            num_haps_inpool[p] = haps.length;
            for (int h = 0; h < haps.length; h++){ // ...for each raw GC haplotype...
                String curr_vc = ""; // the hap_string
                int hap_ct = Integer.parseInt(haps[h].split("\t")[1]);
                for(int k = 0; k < this.num_sites - 1; k++) { // ...for each potential gap...
                    curr_vc += haps[h].charAt(k * 2);
                    // If the linkage between the two variant positions is uncertain
                    // i.e.: gap is present...
                    if (haps[h].charAt(k * 2 + 1) == '?'){
                        gap_counts_inpool[p][k] += hap_ct; // increment the count within-pool and globally
                        gap_counts_all[k] += hap_ct;
                    }
                }
                curr_vc += haps[h].charAt((this.num_sites - 1) * 2);
                // put the haplotype "curr_vc" and its count to the table "hap_tracker"
                if (!hap_tracker.containsKey(curr_vc)) {
                    hap_tracker.put(curr_vc, hap_ct);
                } else {
                    int new_ct = hap_tracker.get(curr_vc) + hap_ct;
                    hap_tracker.put(curr_vc, new_ct);
                }
                tot_ct += hap_ct;
            }
        }
        //Step 2) Check if the potential gaps meet the within- OR between-pool frequency thresholds.
        HashSet<Integer> gap_indexes_set = new HashSet<Integer>();
        for (int k = 0; k < this.num_sites - 1; k++) {
            if ((double) gap_counts_all[k] / num_haps_all >= this.dp.gap_all_pool_cutoff
                && !gap_indexes_set.contains(k)) {

                gap_indexes.add(k);
                gap_indexes_set.add(k);
            }

            for (int p = 0; p < this.num_pools; p++) {
                if ((double) gap_counts_inpool[p][k] / num_haps_inpool[p] >= this.dp.gap_inpool_cutoff
                    && !gap_indexes_set.contains(k)) {

                    gap_indexes.add(k);
                    gap_indexes_set.add(k);
                }
            }
        }
        // Add the last site index as the final "gap" so that all the sites are within gaps.
        gap_indexes.add(this.num_sites - 1);
        // This is not useful for the moment, but keep the data integrity for potential future use.
        gap_indexes_set.add(this.num_sites - 1);

        // Clean the outcome and return.
        int[] gap_indexes_array=new int[gap_indexes.size()];
        for (int i = 0; i < gap_indexes.size(); i++) {
            gap_indexes_array[i] = gap_indexes.get(i);
        }
        Arrays.sort(gap_indexes_array); // sort the indices

        this.global_haps_gc = new String[hap_tracker.size()][this.num_sites];
        this.global_gc_freq = new double[hap_tracker.size()];
        this.num_haps_gc = hap_tracker.size();
        int hap_index = 0;
        for (String curr_vc : hap_tracker.keySet()) {
            String[] tmp = curr_vc.split("");
            for (int l = 0; l < this.num_sites; l++) {
                this.global_haps_gc[hap_index][l] = tmp[l];
            }
            this.global_gc_freq[hap_index] = ((double)hap_tracker.get(curr_vc))/ ((double)tot_ct);
            hap_index++;
        }
        return gap_indexes_array;
    }

    /**
     *  Output the current dividing plan to a file
     *
     *  @param dc_plan_outfile
     */
    public void output_current_DC_plan(String dc_plan_outfile) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(dc_plan_outfile));

            // Output Level I:
            bw.write("Level_I\n");
            System.out.println("Level_I");
            for (int r = 0; r < this.num_regions_level_I; r++) {
                bw.write(this.regions_level_I[r][0] + ":" + this.regions_level_I[r][1] + "\t");
                System.out.print(this.regions_level_I[r][0]
                    + ":"
                    + this.regions_level_I[r][1]
                    + "\t");

            }

            // Output Level II:
            bw.write("\nLevel_II\n");
            System.out.println("\nLevel_II");
            for (int r = 0; r < this.num_regions_level_II; r++) {
                bw.write(this.regions_level_II[r][0] + ":" + this.regions_level_II[r][1] + "\t");
                System.out.print(this.regions_level_II[r][0]
                    + ":"
                    + this.regions_level_II[r][1]
                    + "\t");
            

            }
            bw.write("\nLevel_III\n");
            System.out.println("\nLevel_III");
            for (int r = 0; r < this.num_regions_level_III; r++) {
                bw.write(this.regions_level_III[r][0] + ":" + this.regions_level_III[r][1] + "\t");
                System.out.print(this.regions_level_III[r][0]
                    + ":"
                    + this.regions_level_III[r][1]
                    + "\t");
            }
            bw.write("\nLevel_IV\n");
            System.out.println("\nLevel_IV");
            for (int r = 0; r < this.num_regions_level_IV; r++) {
                bw.write(this.regions_level_IV[r][0] + ":" + this.regions_level_IV[r][1] + "\t");
                System.out.print(this.regions_level_IV[r][0]
                    + ":"
                    + this.regions_level_IV[r][1]
                    + "\t");
            }
            
            bw.write("\nLevel_V\n");
            System.out.println("\nLevel_V");
            for (int r = 0; r < this.num_regions_level_V; r++) {
                bw.write(this.regions_level_V[r][0] + ":" + this.regions_level_V[r][1] + "\t");
                System.out.print(this.regions_level_V[r][0]
                    + ":"
                    + this.regions_level_V[r][1]
                    + "\t");
            }
            bw.write("\nLevel_VI\n");
            System.out.println("\nLevel_VI");
            for (int r = 0; r < this.num_regions_level_VI; r++) {
                bw.write(this.regions_level_VI[r][0] + ":" + this.regions_level_VI[r][1] + "\t");
                System.out.print(this.regions_level_VI[r][0]
                    + ":"
                    + this.regions_level_VI[r][1]
                    + "\t");
            }
            bw.write("\nLevel_VII\n");
            System.out.println("\nLevel_VII");
            for (int r = 0; r < this.num_regions_level_VII; r++) {
                bw.write(this.regions_level_VII[r][0] + ":" + this.regions_level_VII[r][1] + "\t");
                System.out.print(this.regions_level_VII[r][0]
                    + ":"
                    + this.regions_level_VII[r][1]
                    + "\t");
            }
            bw.write("\nLevel_VIII\n");
            System.out.println("\nLevel_VIII");
            for (int r = 0; r < this.num_regions_level_VIII; r++) {
                bw.write(this.regions_level_VIII[r][0] + ":" + this.regions_level_VIII[r][1] + "\t");
                System.out.print(this.regions_level_VIII[r][0]
                    + ":"
                    + this.regions_level_VIII[r][1]
                    + "\t");
            }
            
            
            
            
            
            bw.write("\n");
            
            
            bw.close();
            System.out.println();

        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    
    
    
    double freq_cal (double [][] link_freq, String hap, int start, int end){
    	if (end==start){
    	 return 0.5;
    	}
    	HashMap<String, Integer> hap2index= new HashMap<String, Integer>();
    	hap2index.put("00",0);
    	hap2index.put("01",1);
    	hap2index.put("10",2);
    	hap2index.put("11",3);
    	double min_freq =1.0;
    	double cal_freq=0.0;
    	for (int i=0; i<(end-start-1); i++ ){
    		String sub_hap= hap.substring(i, i+2);  
    		double freq= link_freq[start+i][hap2index.get(sub_hap)];
    		if (freq> 0.001) {
    			if (freq< min_freq){
    				min_freq= freq;
    			}
    			cal_freq+= freq;
    		}else{
    			return 0.0;
    		}
    	}
    	return min_freq;
//    	return cal_freq;
    }

    /**
     *  Analyze a set of regions based on the dividing plan.
     *
     *  This function can be used for either level I or level II, or a subset of them
     *
     *  @param regions
     *  @param parameter_file
     *  @param dir_prefix
     *  @param level_index
     *  @return
     *  @throws Exception
     */
    
    public HapConfig[] regional_AEM(
    	String[] pool_IDs,
    	String[] vef_files,
        int[][] regions,
        String parameter_file,
        String dir_prefix,
        String aem_fail_lasso_path,
        int level_index) throws Exception {

        HapConfig[] region_haps = new HapConfig[regions.length];
        for (int r_index = 0; r_index < regions.length; r_index++) {
            System.out.print("AEM for level " + level_index + " region " + r_index + "... ");
            HapConfig hap_config = generate_hapconfig_2n(regions, r_index, pool_IDs);
            
            RegionEMSolver hap_solver = new RegionEMSolver(hap_config, parameter_file,0 ,r_index,1 );
            if (hap_solver.failure) {
                System.out.println("AEM failed to converge. Initiating regional LASSO... ");

                // If AEM fails, at least some of the frequencies will be NaN. In that case, use
                // sub-optimal GC results. TODO [Quan] 
                region_haps[r_index] = generate_hapconfig_gc_regional_lasso(pool_IDs, vef_files,
                    regions[r_index], level_index, r_index, aem_fail_lasso_path);
            } else {
                region_haps[r_index] = hap_solver.final_Haps;
            }
            region_haps[r_index].recode_HapIDs_to_base16();
            region_haps[r_index].write_global_file_string(dir_prefix
                + "_level_"
                + level_index
                + "_region_"
                + r_index
                + ".inter_freq_haps.txt");

            System.out.print("Done. AEM ");
            if (hap_solver.failure) {
                System.out.print("failed to converge. ");
            } else {
                System.out.print("successfully converged. ");
            }

            System.out.println(region_haps[r_index].num_global_hap
                + " regional haplotypes have been generated.");
        }

        return region_haps;
    }
    
    /**
     * @author Chen Cao 2019-08
     * Using the frequency path to  Evaluate whether the haplotypes ( 2^N )exist 
     * And Estimate the frequency of the haplotypes.
     * The output is used as the initial matrix for AEM 
     * 
     */  

    
    public HapConfig[] regional_AEM(
        	String[] pool_IDs,
        	String[] vef_files,
            int[][] regions,
            String parameter_file,
            String dir_prefix,
            String aem_fail_lasso_path,
            int level_index, 
            HashMap<String, Integer> name2index) throws Exception {
    		
            HapConfig[] region_haps = new HapConfig[regions.length];
            for (int r_index = 0; r_index < regions.length; r_index++) {
                System.out.print("AEM for level " + level_index + " region " + r_index + "... ");
                HapConfig hap_config = generate_hapconfig_gc_exhaustive(regions, r_index,
                		pool_IDs, name2index);
                
                RegionEMSolver hap_solver = new RegionEMSolver(hap_config, parameter_file, 1, r_index, 1 );
                if (hap_solver.failure) {
                	System.exit(0);
                	double  final_diff =1.0;
                	
                	HapConfig hap_config_all = generate_hapconfig_2n(regions, r_index,
                    		pool_IDs);
                	RegionEMSolver hap_solver_all= new 
                			RegionEMSolver(hap_config_all, parameter_file, 1, r_index, 1 );
                	
                	if (hap_solver_all.failure) {
                		if (hap_solver_all.aem_min_diff< final_diff) {
                			final_diff= hap_solver_all.aem_min_diff ;
                		}
                		RegionEMSolver hap_solver_fail = new RegionEMSolver(hap_config_all, parameter_file, -1, 
                				-1, final_diff+this.aem_fail_dist_cutoff);
                		region_haps[r_index] = hap_solver_fail.final_Haps;
                	}else {
                		region_haps[r_index] = hap_solver_all.final_Haps;
                	}	
                } else {
                    region_haps[r_index] = hap_solver.final_Haps;
                }
                region_haps[r_index].recode_HapIDs_to_base16();
                region_haps[r_index].write_global_file_string(dir_prefix
                    + "_level_"
                    + level_index
                    + "_region_"
                    + r_index
                    + ".inter_freq_haps.txt");

                System.out.print("Done. AEM ");
                if (hap_solver.failure) {
                    System.out.print("failed to converge. ");
                } else {
                    System.out.print("successfully converged. ");
                }

                System.out.println(region_haps[r_index].num_global_hap
                    + " regional haplotypes have been generated.");
            }

            return region_haps;
    }
    
    public HapConfig[] level_III_regional_AEM(
        	String[] pool_IDs,
        	String[] vef_files,
            int[][] regions,
            String parameter_file,
            String dir_prefix,
            String aem_fail_lasso_path,
            int level_index, 
            HashMap<String, Integer> name2index,
            int mismatch_tolerance) throws Exception {
    		
            HapConfig[] region_haps = new HapConfig[regions.length];
            this.level_III_IV_region_mismatch_tolerance= mismatch_tolerance;
            for (int r_index = 0; r_index < regions.length; r_index++) {
//            for (int r_index = 2; r_index < regions.length; r_index++) { 	
                System.out.print("AEM for level " + level_index + " region " + r_index + "... ");
                
                HapConfig hap_config = generate_hapconfig_level_III(regions, r_index,
                		pool_IDs);
                
                
                RegionEMSolver hap_solver = new RegionEMSolver(hap_config, parameter_file, 3,r_index, 1);
                
                
                if (hap_solver.failure) {
                	int final_index= -1;
                	double  final_diff =1.0;
                	boolean all_fail= true;
                	
                	int index=-1;
                	while (index<this.aem_fail_max_iter) {
	                	HapConfig hap_config_fail = generate_hapconfig_level_III_fail(regions, r_index,
	                        		pool_IDs, index);
	                	RegionEMSolver hap_solver_fail = new RegionEMSolver(hap_config_fail, parameter_file, 3, r_index,1);
	                	if  (!hap_solver_fail.failure) {
	                		all_fail =false;
	                		region_haps[r_index] = hap_solver_fail.final_Haps;
	                		
	                		break;
	                	}else {
	                		if (hap_solver_fail.aem_min_diff < final_diff) {
	                			final_diff= hap_solver_fail.aem_min_diff ;
	                			final_index= index;
	                		}
	                	}
                		index++;
                	}
                	
                	if ( all_fail) {
                		HapConfig hap_config_fail = generate_hapconfig_level_III_fail(regions, r_index,
                        		pool_IDs, final_index);
                		RegionEMSolver hap_solver_fail = new RegionEMSolver(hap_config_fail, parameter_file, -1, 
                				-1, final_diff+this.aem_fail_dist_cutoff);
                		region_haps[r_index] = hap_solver_fail.final_Haps;
                	}
                	
                	
//                    System.out.println("AEM failed to converge. Initiating regional LASSO... ");

                    // If AEM fails, at least some of the frequencies will be NaN. In that case, use
                    // sub-optimal GC results. TODO [Quan] 
//                    region_haps[r_index] = generate_hapconfig_gc_regional_lasso(pool_IDs, vef_files,
//                        regions[r_index], level_index, r_index, aem_fail_lasso_path);
                    
                } else {
                    region_haps[r_index] = hap_solver.final_Haps;
                }
                region_haps[r_index].recode_HapIDs_to_base16();
                region_haps[r_index].write_global_file_string(dir_prefix
                    + "_level_"
                    + level_index
                    + "_region_"
                    + r_index
                    + ".inter_freq_haps.txt");

                System.out.print("Done. AEM ");
//                if (hap_solver.failure) {
//                    System.out.print("failed to converge. ");
//                } else {
//                    System.out.print("successfully converged. ");
//                }

                System.out.println(region_haps[r_index].num_global_hap
                    + " regional haplotypes have been generated.");
            }

            return region_haps;
    }
    
    public HapConfig[] level_IV_regional_AEM(
        	String[] pool_IDs,
        	String[] vef_files,
            int[][] regions,
            String parameter_file,
            String dir_prefix,
            String aem_fail_lasso_path,
            int level_index, 
            HashMap<String, Integer> name2index,
            int mismatch_tolerance) throws Exception {
    		
    		this.level_III_IV_region_mismatch_tolerance= mismatch_tolerance;
            HapConfig[] region_haps = new HapConfig[regions.length];
            for (int r_index = 0; r_index < regions.length; r_index++) {
                System.out.print("AEM for level " + level_index + " region " + r_index + "... ");
                
                HapConfig hap_config = generate_hapconfig_level_IV(regions, r_index,
                		pool_IDs);
                
//        		String file_path= "/home/chencao/Desktop/cc.txt";
//        		FileWriter mydata = new FileWriter(file_path,false);
//                PrintWriter pw = new PrintWriter(mydata);
//                String pw_line ="";
//               
//                pw.flush();
//                pw.close();
                
                RegionEMSolver hap_solver = new RegionEMSolver(hap_config, parameter_file,4 ,r_index,1 );
                
                if (hap_solver.failure) {
                	int final_index= -1;
                	double  final_diff =1.0;
                	boolean all_fail= true;
                	
                	int index=-1;
                	while (index<this.aem_fail_max_iter) {
	                	HapConfig hap_config_fail = generate_hapconfig_level_IV_fail(regions, r_index,
	                        		pool_IDs, index);
	                	RegionEMSolver hap_solver_fail = new RegionEMSolver(hap_config_fail, parameter_file, 4,r_index, 1);
	                	if  (!hap_solver_fail.failure) {
	                		all_fail =false;
	                		region_haps[r_index] = hap_solver_fail.final_Haps;
	                		break;
	                	}else {
	                		if (hap_solver_fail.aem_min_diff < final_diff) {
	                			final_diff= hap_solver_fail.aem_min_diff ;
	                			final_index= index;
	                		}
	                	}
                		index++;
                	}
                	if ( all_fail) {
                		HapConfig hap_config_fail = generate_hapconfig_level_IV_fail(regions, r_index,
                        		pool_IDs, final_index);
                		RegionEMSolver hap_solver_fail = new RegionEMSolver(hap_config_fail, parameter_file, -1, 
                				-1, final_diff+this.aem_fail_dist_cutoff);
                		region_haps[r_index] = hap_solver_fail.final_Haps;
                	}
                    
                } else {
                    region_haps[r_index] = hap_solver.final_Haps;
                }
                region_haps[r_index].recode_HapIDs_to_base16();
                region_haps[r_index].write_global_file_string(dir_prefix
                    + "_level_"
                    + level_index
                    + "_region_"
                    + r_index
                    + ".inter_freq_haps.txt");

                System.out.print("Done. AEM ");
//                if (hap_solver.failure) {
//                    System.out.print("failed to converge. ");
//                } else {
//                    System.out.print("successfully converged. ");
//                }

                System.out.println(region_haps[r_index].num_global_hap
                    + " regional haplotypes have been generated.");
            }

            return region_haps;
    }
    
    public HapConfig[] level_V_regional_AEM(
        	String[] pool_IDs,
        	String[] vef_files,
            int[][] regions,
            String parameter_file,
            String dir_prefix,
            String aem_fail_lasso_path,
            int level_index, 
            HashMap<String, Integer> name2index,
            int mismatch_tolerance) throws Exception {
    		this.level_V_VI_region_mismatch_tolerance = mismatch_tolerance;
            HapConfig[] region_haps = new HapConfig[regions.length];
            for (int r_index = 0; r_index < regions.length; r_index++) {
                System.out.print("AEM for level " + level_index + " region " + r_index + "... ");
                HapConfig hap_config = generate_hapconfig_level_V(regions, r_index,
                		pool_IDs);
                RegionEMSolver hap_solver = new RegionEMSolver(hap_config, parameter_file, 5,r_index, 1 );
                if (hap_solver.failure) {
                	
                	int final_index= -1;
                	double  final_diff =1.0;
                	boolean all_fail= true;
                	
                	int index=-2;
                	while (index<this.aem_fail_max_iter) {
	                	HapConfig hap_config_fail = generate_hapconfig_level_V_fail(regions, r_index,
	                        		pool_IDs, index);
	                	RegionEMSolver hap_solver_fail = new RegionEMSolver(hap_config_fail, parameter_file, 5,r_index, 1);
	                	if  (!hap_solver_fail.failure) {
	                		all_fail =false;
	                		region_haps[r_index] = hap_solver_fail.final_Haps;
	                		
	                		break;
	                	}else {
	                		if (hap_solver_fail.aem_min_diff < final_diff) {
	                			final_diff= hap_solver_fail.aem_min_diff ;
	                			final_index= index;
	                		}
	                	}
                		index++;
                	}
                	System.out.println(final_diff);
                	if ( all_fail) {
                		HapConfig hap_config_fail = generate_hapconfig_level_V_fail(regions, r_index,
                        		pool_IDs, final_index);
                		RegionEMSolver hap_solver_fail = new RegionEMSolver(hap_config_fail, parameter_file, -1, 
                				-1, final_diff+this.aem_fail_dist_cutoff);
                		region_haps[r_index] = hap_solver_fail.final_Haps;
                	}
                    
                } else {
                    region_haps[r_index] = hap_solver.final_Haps;
                }
                region_haps[r_index].recode_HapIDs_to_base16();
                region_haps[r_index].write_global_file_string(dir_prefix
                    + "_level_"
                    + level_index
                    + "_region_"
                    + r_index
                    + ".inter_freq_haps.txt");

                System.out.print("Done. AEM ");
//                if (hap_solver.failure) {
//                    System.out.print("failed to converge. ");
//                } else {
//                    System.out.print("successfully converged. ");
//                }

                System.out.println(region_haps[r_index].num_global_hap
                    + " regional haplotypes have been generated.");
            }
            return region_haps;
    }
    
    public HapConfig[] level_VI_regional_AEM(
        	String[] pool_IDs,
        	String[] vef_files,
            int[][] regions,
            String parameter_file,
            String dir_prefix,
            String aem_fail_lasso_path,
            int level_index, 
            HashMap<String, Integer> name2index,
            int mismatch_tolerance) throws Exception {
    		this.level_V_VI_region_mismatch_tolerance = mismatch_tolerance;
            HapConfig[] region_haps = new HapConfig[regions.length];
            for (int r_index = 0; r_index < regions.length; r_index++) {
                System.out.print("AEM for level " + level_index + " region " + r_index + "... ");
                HapConfig hap_config = generate_hapconfig_level_VI(regions, r_index,
                		pool_IDs);
                RegionEMSolver hap_solver = new RegionEMSolver(hap_config, parameter_file,6,r_index,1 );
                if (hap_solver.failure) {
                	int final_index= -1;
                	double  final_diff =1.0;
                	boolean all_fail= true;
                	
                	int index=-2;
                	while (index<this.aem_fail_max_iter) {
	                	HapConfig hap_config_fail = generate_hapconfig_level_VI_fail(regions, r_index,
	                        		pool_IDs, index);
	                	RegionEMSolver hap_solver_fail = new RegionEMSolver(hap_config_fail, parameter_file, 6, r_index, 1);
	                	if  (!hap_solver_fail.failure) {
	                		all_fail =false;
	                		region_haps[r_index] = hap_solver_fail.final_Haps;
	                		
	                		break;
	                	}else {
	                		if (hap_solver_fail.aem_min_diff < final_diff) {
	                			final_diff= hap_solver_fail.aem_min_diff ;
	                			final_index= index;
	                		}
	                	}
                		index++;
                	}
                	if ( all_fail) {
                		HapConfig hap_config_fail = generate_hapconfig_level_VI_fail(regions, r_index,
                        		pool_IDs, final_index);
                		RegionEMSolver hap_solver_fail = new RegionEMSolver(hap_config_fail, parameter_file, -1, 
                				-1 ,final_diff+this.aem_fail_dist_cutoff);
                		region_haps[r_index] = hap_solver_fail.final_Haps;
                	}
                    
                } else {
                    region_haps[r_index] = hap_solver.final_Haps;
                }
                region_haps[r_index].recode_HapIDs_to_base16();
                region_haps[r_index].write_global_file_string(dir_prefix
                    + "_level_"
                    + level_index
                    + "_region_"
                    + r_index
                    + ".inter_freq_haps.txt");

                System.out.print("Done. AEM ");
                System.out.println(region_haps[r_index].num_global_hap
                    + " regional haplotypes have been generated.");
            }

            return region_haps;
    }
    
    
    public HapConfig[] level_VII_regional_AEM(
        	String[] pool_IDs,
        	String[] vef_files,
            int[][] regions,
            String parameter_file,
            String dir_prefix,
            String aem_fail_lasso_path,
            int level_index, 
            HashMap<String, Integer> name2index,
            int mismatch_tolerance) throws Exception {
    		this.level_VII_VIII_region_mismatch_tolerance = mismatch_tolerance;
            HapConfig[] region_haps = new HapConfig[regions.length];
            for (int r_index = 0; r_index < regions.length; r_index++) {
                System.out.print("AEM for level " + level_index + " region " + r_index + "... ");
                HapConfig hap_config = this.generate_hapconfig_level_VII (regions, r_index,
                		pool_IDs);
                RegionEMSolver hap_solver = new RegionEMSolver(hap_config, parameter_file, 7,r_index, 1 );
                if (hap_solver.failure) {
                	
                	int final_index= -1;
                	double  final_diff =1.0;
                	boolean all_fail= true;
                	
                	int index=-2;
                	while (index<(this.aem_fail_max_iter+1)) {
	                	HapConfig hap_config_fail = generate_hapconfig_level_VII_fail(regions, r_index,
	                        		pool_IDs, index);
	                	RegionEMSolver hap_solver_fail = new 
	                			RegionEMSolver(hap_config_fail, parameter_file, 7,r_index, 1);
	                	if  (!hap_solver_fail.failure) {
	                		all_fail =false;
	                		region_haps[r_index] = hap_solver_fail.final_Haps;
	                		
	                		break;
	                	}else {
	                		if (hap_solver_fail.aem_min_diff < final_diff) {
	                			final_diff= hap_solver_fail.aem_min_diff ;
	                			final_index= index;
	                		}
	                	}
                		index++;
                	}
                	System.out.println(final_diff);
                	if ( all_fail) {
                		HapConfig hap_config_fail = generate_hapconfig_level_VII_fail(regions, r_index,
                        		pool_IDs, final_index);
                		RegionEMSolver hap_solver_fail = new RegionEMSolver(hap_config_fail, parameter_file, -1, 
                				-1, final_diff+this.aem_fail_dist_cutoff);
                		region_haps[r_index] = hap_solver_fail.final_Haps;
                	}
                    
                } else {
                    region_haps[r_index] = hap_solver.final_Haps;
                }
                region_haps[r_index].recode_HapIDs_to_base16();
                region_haps[r_index].write_global_file_string(dir_prefix
                    + "_level_"
                    + level_index
                    + "_region_"
                    + r_index
                    + ".inter_freq_haps.txt");

                System.out.print("Done. AEM ");
//                if (hap_solver.failure) {
//                    System.out.print("failed to converge. ");
//                } else {
//                    System.out.print("successfully converged. ");
//                }

                System.out.println(region_haps[r_index].num_global_hap
                    + " regional haplotypes have been generated.");
            }
            return region_haps;
    }
    
    
    public HapConfig[] level_VIII_regional_AEM(
        	String[] pool_IDs,
        	String[] vef_files,
            int[][] regions,
            String parameter_file,
            String dir_prefix,
            String aem_fail_lasso_path,
            int level_index, 
            HashMap<String, Integer> name2index,
            int mismatch_tolerance) throws Exception {
    		this.level_VII_VIII_region_mismatch_tolerance = mismatch_tolerance;
            HapConfig[] region_haps = new HapConfig[regions.length];
            for (int r_index = 0; r_index < regions.length; r_index++) {
                System.out.print("AEM for level " + level_index + " region " + r_index + "... ");
                HapConfig hap_config = generate_hapconfig_level_VIII(regions, r_index,
                		pool_IDs);
                RegionEMSolver hap_solver = new RegionEMSolver(hap_config, parameter_file,8,r_index,1 );
                if (hap_solver.failure) {
                	int final_index= -1;
                	double  final_diff =1.0;
                	boolean all_fail= true;
                	
                	int index=-2;
                	while (index<(this.aem_fail_max_iter+1)) {
	                	HapConfig hap_config_fail = generate_hapconfig_level_VIII_fail(regions, r_index,
	                        		pool_IDs, index);
	                	
		                	RegionEMSolver hap_solver_fail = new RegionEMSolver(hap_config_fail, parameter_file, 8, r_index, 1);
		                	if  (!hap_solver_fail.failure) {
		                		all_fail =false;
		                		region_haps[r_index] = hap_solver_fail.final_Haps;
		                		
		                		break;
		                	}else {
		                		if (hap_solver_fail.aem_min_diff < final_diff) {
		                			final_diff= hap_solver_fail.aem_min_diff ;
		                			final_index= index;
		                		}
		                	}
                		index++;
                	}
                	if ( all_fail) {
                		HapConfig hap_config_fail = generate_hapconfig_level_VIII_fail(regions, r_index,
                        		pool_IDs, final_index);
                		RegionEMSolver hap_solver_fail = new RegionEMSolver(hap_config_fail, parameter_file, -1, 
                				-1 ,final_diff+ this.aem_fail_dist_cutoff);
                		region_haps[r_index] = hap_solver_fail.final_Haps;
                	}
                    
                } else {
                    region_haps[r_index] = hap_solver.final_Haps;
                }
                region_haps[r_index].recode_HapIDs_to_base16();
                region_haps[r_index].write_global_file_string(dir_prefix
                    + "_level_"
                    + level_index
                    + "_region_"
                    + r_index
                    + ".inter_freq_haps.txt");

                System.out.print("Done. AEM ");
                System.out.println(region_haps[r_index].num_global_hap
                    + " regional haplotypes have been generated.");
            }

            return region_haps;
    }
    
    
    public boolean Is2Power(int x, int n) {
    	if (x==0) {
    		return true;
    	}
    	for (int i = 0; i < n; i++) {
    		if (x== Math.pow(2, i)){
    			return true;
    		}
    	}
    	return false;
    }
    
    public boolean Is2PowerMinus(int x, int n) {

    	for (int i = 0; i < n; i++) {
    		if (x==  (Math.pow(2, n)-1-  Math.pow(2, i)) ){
    			return true;
    		}
    	}
    	return false;
    }
    
    public int num_of_alternate (String [] hap ) throws IOException {
    	int count= 0;
    	for (int i =0; i< hap.length;i++) {
    		if (hap[i].equals("1")) {
    			 count++;
    		}
    	}
    	return count;
    }
    
    public String [] strmatch (String x, String y, String z, 
    		int num, int num_mismatch_cutoff)  throws IOException {
		
    	ArrayList<Integer >  pos_arr = new ArrayList<Integer>();
//    	System.out.println( x );
//    	System.out.println( y );
//    	System.out.println( z );
//    	System.out.println( num );
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
//    	System.out.println( overlap_y );
//    	System.out.println( overlap_yz );
    	for (int i = 0; i < overlap_y.length(); i++) {
			if (overlap_y.charAt(i) != overlap_yz.charAt(i)) {
				distance++;
				pos_arr.add(x.length()+ i);
			}
		}
    	if (distance > num_mismatch_cutoff) {
    		String[] com_str2 = new String [ 0];
    		return com_str2;
    	}
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
        	// the length of vc_str may not reach num_site_regional; so put zeros in.
        	for (int locus = 0; locus < distance - l ; locus++) {
        		vc_str = "0"+ vc_str;
            }
//        	System.out.println( vc_str);
        	sub_str[h] = vc_str;
    	}
    	
    	for (int h = 0; h < haps_2n; h++) {
    		com_str[h]= link_str(com_str[h], sub_str[h], pos_arr );
//    		System.out.println( com_str[h] );
    	}
    	
    	return com_str;
    }
    
    public  String link_str(String x, String y,ArrayList<Integer >  pos_arr  )
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
    
    public String diagonal_str (int index, int n) 
    		throws IOException {
    	String tmp = "";
    	for (int i=0;i< n;i++) {
    		if (i==index) {
    			tmp=tmp+"1";
    		}else {
    			tmp=tmp+"0";
    		}
    	}
    	return tmp;
    }
    
    
    public String diagonal_0_str (int index, int n) 
    		throws IOException {
    	String tmp = "";
    	for (int i=0;i< n;i++) {
    		if (i==index) {
    			tmp=tmp+"0";
    		}else {
    			tmp=tmp+"1";
    		}
    	}
    	return tmp;
    }

    /**
     * @author Chen Cao 2019-08
     * Construct the frequency path for four possible ways: 0/0 0/1 1/0 1/1
     * Using the frequency path to  Evaluate whether the haplotypes ( 2^N )exist 
     * And Estimate the frequency of the haplotypes.
     * Using exhaustive to find out all possible haplotypes.
     * 
     */  
    
    public HapConfig generate_hapconfig_gc_exhaustive(
    	int[][] the_region, int region_index, String[] pool_IDs, 
    		HashMap<String, Integer> name2index_dict) throws IOException {
    		
    		
            int region_start = the_region[region_index][0];
            int region_end = the_region[region_index][1];
            int num_site_regional = region_end-region_start + 1;
            double[][] inpool_site_freqs = new double[num_site_regional][];
            LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];
            for (int l = 0; l < num_site_regional; l++) {
                locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
                inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
            }
            boolean Using1asDiagonal=true;
            double inpool_site_total_freq=0.0;
            double inpool_site_total_count=0.0;
            for (int l = 0; l < num_site_regional; l++) {
            	for (int i=0; i<inpool_site_freqs[l].length; i++ ) {
            		inpool_site_total_freq+= inpool_site_freqs[l][i];
            		inpool_site_total_count+=1.0;
            	}
            }
            if ((inpool_site_total_freq / inpool_site_total_count) <0.5) {
            	Using1asDiagonal=false;
            }
            System.out.println((inpool_site_total_freq) /inpool_site_total_count);
            
            int haps_2n = (int) Math.pow(2, num_site_regional);
            String[][] global_haps_string_tmp  = new String[haps_2n][num_site_regional];
            String[] hap_IDs_tmp = new String[haps_2n];
            double[] global_haps_freq_tmp  = new double [haps_2n];
            double total_freq =0;
            for (int h = 0; h < haps_2n; h++) {
            	String curr_ID = "";
            	String vc_str = Integer.toBinaryString(h);
            	String[] vc_arr  = vc_str.split("");
            	// the length of vc_str may not reach num_site_regional; so put zeros in.
            	for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                    global_haps_string_tmp[h][locus] = "0";
                    curr_ID += "0";
                }
            	for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
            		global_haps_string_tmp[h][locus] = vc_arr[locus - num_site_regional 
            		                                          + vc_arr.length];
            		curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
            	}
            	hap_IDs_tmp[h] = curr_ID;
            	for (int p=0; p< this.num_pools;p++) {
            		double freq = freq_cal(this.loci_link_freq[p],curr_ID, region_start,region_end );
            		global_haps_freq_tmp[h]+= freq ;
            		total_freq+= freq;
            	}
            }
            
            int count=0;
            total_freq=total_freq+ (double)num_site_regional*0.001;
            
            
            boolean [] AsDiagonal = new boolean[global_haps_freq_tmp.length];
            for (int i=0; i < global_haps_freq_tmp.length;i++) {
            	AsDiagonal[i] =false;
            	boolean i_ok=false;
            	if ((global_haps_freq_tmp[i]/ total_freq) > (0.1/ (double) haps_2n)) {
            		i_ok = true;
            	}
//            	if ((Using1asDiagonal) && Is2Power(i, num_site_regional)) {
//            		i_ok = true;
//            		AsDiagonal[i] =true;
//        		}
//        		if ((!Using1asDiagonal) && Is2PowerMinus(i, num_site_regional)) {
//        			i_ok = true;
//        			AsDiagonal[i] =true;
//        		}
        		if (i_ok) {
        			count++;
        		}
            }
            
            double power2_min= 0.001; 
            int fix_count= count;
            String[][] global_haps_string = new String[count][num_site_regional];
            String[] hap_IDs = new String[count];
            double[] global_haps_freq = new double[count];
            count=0;
            for (int i=0; i < global_haps_freq_tmp.length;i++) {
            	if (((global_haps_freq_tmp[i]/ total_freq) > (0.1 / (double) haps_2n))  || 
            			AsDiagonal[i] )	{
	            	for (int j=0; j < global_haps_string[count].length;j++) {
	            		global_haps_string[count][j]= global_haps_string_tmp[i][j];
	           		}
	           		hap_IDs[count]= hap_IDs_tmp[i];
	           		global_haps_freq[count]= global_haps_freq_tmp[i]/ total_freq;
	           		if (AsDiagonal[i]) {
	           			global_haps_freq[count]= (global_haps_freq_tmp[i]+power2_min )/ total_freq;
	           		}

	           		count++;
            	
            	}
            }

            
            return new HapConfig(
                global_haps_string,
                global_haps_freq,
                null,
                inpool_site_freqs,
                locusInfo,
                this.num_pools,
                hap_IDs,
                pool_IDs,
                this.dp.est_ind_pool,
                Using1asDiagonal);    	
    }
    
    /**
     *	generate the HapConfig with 2^n haplotypes (n is the number of sites).
     *
     *  @param the_region
     *  @param region_index
     *  @param pool_IDs
     *  @return
     *  @throws IOException
     */	
    
    public HapConfig generate_hapconfig_2n(
        int[][] the_region, int region_index, String[] pool_IDs) throws IOException {

        int region_start = the_region[region_index][0];
        int region_end = the_region[region_index][1];
        int num_site_regional = region_end-region_start + 1;
        double[][] inpool_site_freqs = new double[num_site_regional][];
        LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];

        for (int l = 0; l < num_site_regional; l++) {
            locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
            inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
        }
        int haps_2n = (int) Math.pow(2, num_site_regional);
        String[][] global_haps_string = new String[haps_2n][num_site_regional];
        String[] hap_IDs = new String[haps_2n];
        for (int h = 0; h < haps_2n; h++) {
            String curr_ID = "";
            String vc_str = Integer.toBinaryString(h);
            String[] vc_arr  = vc_str.split("");
            // the length of vc_str may not reach num_site_regional; so put zeros in.
            for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                global_haps_string[h][locus] = "0";
                curr_ID += "0";
            }
            for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
                global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
                curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
            }
            hap_IDs[h] = curr_ID;
        }

        double[] global_haps_freq = new double[haps_2n];
        Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
        return new HapConfig(
            global_haps_string,
            global_haps_freq,
            null,
            inpool_site_freqs,
            locusInfo,
            this.num_pools,
            hap_IDs,
            pool_IDs,
            this.dp.est_ind_pool);
    }
    

    
    public HapConfig generate_hapconfig_level_III(
            int[][] the_region, int region_index, String[] pool_IDs) throws IOException {

            int region_start = the_region[region_index][0];
            int region_end = the_region[region_index][1];
            int num_site_regional = region_end-region_start + 1;
            double[][] inpool_site_freqs = new double[num_site_regional][];
            LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];

            for (int l = 0; l < num_site_regional; l++) {
                locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
                inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
            }
            HashSet<String > haps_set =new HashSet<String >();
            ArrayList<String >  linked_haps= new ArrayList<String >();
            
//            if (num_site_regional==28) {
//            	linked_haps.add("0110100000010100000001000001");
//            	linked_haps.add("1000001111000001110100101100");
//            	linked_haps.add("0000001000011100000001000001");
//            	linked_haps.add("1000001111000001110100101100");
//
//            }
            
            
            if ((this.num_regions_level_I%2==1) && (region_index==this.num_regions_level_III-1)){
            	for (int i =0;i< this.level_I_config[ this.level_I_config.length-1].num_global_hap;
            			i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_I_config[ this.level_I_config.length-1].
            				global_haps_string[i].length;j++) {
            			tmp=tmp+ this.level_I_config[ this.level_I_config.length-1].
                				global_haps_string[i][j];
            		}
            		linked_haps.add(tmp);
            	}
            }else {
            	ArrayList<String >  String_X= new ArrayList<String >();
            	ArrayList<String >  String_Y= new ArrayList<String >();
            	ArrayList<String >  String_Z= new ArrayList<String >();
            	int overlap= 1+ this.regions_level_I[2*region_index ][1]- 
            			this.regions_level_II[2*region_index ][0];
//            	
            	
            	for (int i =0;i< this.level_I_config[2*region_index].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_I_config[2*region_index].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_I_config[2*region_index].global_haps_string[i][j];
            		}
            		String_X.add(tmp);
            	}
            	for (int i =0;i< this.level_I_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_I_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_I_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_Y.add(tmp);
            	}
            	
            	for (int i =0;i< this.level_II_config[2*region_index].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_II_config[2*region_index].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_II_config[2*region_index].global_haps_string[i][j];
            		}
            		String_Z.add(tmp);
            	}
//            	
            	for (int i =0;i< String_X.size();i++) {
            		for (int j =0;j< String_Y.size();j++) {
            			for (int k =0;k< String_Z.size();k++) {
            				String[] comb_str = strmatch( String_X.get(i),
            						String_Y.get(j), String_Z.get(k),overlap,
            						this.level_III_IV_region_mismatch_tolerance);
            				
            				
            				if   ((comb_str.length>0) &&
            						(!haps_set.contains(String_X.get(i)+  String_Y.get(j)))) {
            					linked_haps.add(String_X.get(i)+  String_Y.get(j));
        						haps_set.add(String_X.get(i)+  String_Y.get(j));
        						break;
            				}
            				
            				for (int h=0;h<comb_str.length;h++ ) {
            					if (!haps_set.contains(comb_str[h])) {
            						linked_haps.add(comb_str[h]);
            						haps_set.add(comb_str[h]);
            					}
            				}
            			}
            		}
            	}
            }
            
//            for (int i=0;i<linked_haps.size();i++ ) {
//            	System.out.println(linked_haps.get(i).length()); 
//            }
            
            
            
//          
            boolean Using1asDiagonal=true;
            double inpool_site_total_freq=0.0;
            double inpool_site_total_count=0.0;
            for (int l = 0; l < num_site_regional; l++) {
            	for (int i=0; i<inpool_site_freqs[l].length; i++ ) {
            		inpool_site_total_freq+= inpool_site_freqs[l][i];
            		inpool_site_total_count+=1.0;
            	}
            }
            if ((inpool_site_total_freq / inpool_site_total_count) <0.5) {
            	Using1asDiagonal=false;
            }
            System.out.println((inpool_site_total_freq) /inpool_site_total_count);
            
            
            
          //Chen: Construct a diagonal matrix
            
//            for (int i=0;i< num_site_regional;i++) {
//	            	if (Using1asDiagonal) {
//		            	String diagonal_hap =  diagonal_str(i, num_site_regional);
//		            	if (!haps_set.contains(diagonal_hap )) {
//		            		linked_haps.add(diagonal_hap);
//		            	}
//	            	}
//	            	if (!Using1asDiagonal) {
//		            	String diagonal_hap =  diagonal_0_str(i, num_site_regional);
//		            	if (!haps_set.contains(diagonal_hap )) {
//		            		linked_haps.add(diagonal_hap);
//		            	}
//	            	}
//            }
//            Collections.shuffle(linked_haps);
            
            int haps_2n = linked_haps.size();
            String[][] global_haps_string = new String[haps_2n][num_site_regional];
            String[] hap_IDs = new String[haps_2n];
            for (int h = 0; h < haps_2n; h++) {
                String curr_ID = "";
                String vc_str = linked_haps.get(h);
                String[] vc_arr  = vc_str.split("");
                // the length of vc_str may not reach num_site_regional; so put zeros in.
                for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                    global_haps_string[h][locus] = "0";
                    curr_ID += "0";
                }
                for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
                    global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
                    curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
                }
                hap_IDs[h] = curr_ID;
            }
//            for (int i=0;i<global_haps_string.length;i++ ) {
//            	for (int j=0;j <global_haps_string[i].length;j++ ) {
//            		System.out.print(global_haps_string[i][j]); 
//            	}
//            	System.out.println();
//            }
            //Chen: Construct a diagonal matrix
            
            
            
            
            
            double  diagonal_freq = 0.0002;
            double[] global_haps_freq = new double[haps_2n];
            Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
            
            for (int i=0; i<haps_2n;i++ ) {
            	if ((Using1asDiagonal  ) && (num_of_alternate(global_haps_string[i])==1  )) {
            			global_haps_freq[i]= diagonal_freq;
            	}else if  ((!Using1asDiagonal  ) && (num_of_alternate(global_haps_string[i])
            			==( global_haps_string[i].length-1)  )) {
            		global_haps_freq[i]= diagonal_freq;
            	}
            }

            
            return new HapConfig(
                global_haps_string,
                global_haps_freq,
                null,
                inpool_site_freqs,
                locusInfo,
                this.num_pools,
                hap_IDs,
                pool_IDs,
                this.dp.est_ind_pool,
                Using1asDiagonal);
        }
    
    
    public HapConfig generate_hapconfig_level_III_fail(
            int[][] the_region, int region_index, String[] pool_IDs, int flag_index) throws IOException {

            int region_start = the_region[region_index][0];
            int region_end = the_region[region_index][1];
            int num_site_regional = region_end-region_start + 1;
            double[][] inpool_site_freqs = new double[num_site_regional][];
            LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];

            for (int l = 0; l < num_site_regional; l++) {
                locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
                inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
            }
            HashSet<String > haps_set =new HashSet<String >();
            ArrayList<String >  linked_haps= new ArrayList<String >();
            
//            if (num_site_regional==28) {
//            	linked_haps.add("0110100000010100000001000001");
//            	linked_haps.add("1000001111000001110100101100");
//            	linked_haps.add("0000001000011100000001000001");
//            	linked_haps.add("1000001111000001110100101100");
//
//            }
            
            
            if ((this.num_regions_level_I%2==1) && (region_index==this.num_regions_level_III-1)){
            	for (int i =0;i< this.level_I_config[ this.level_I_config.length-1].num_global_hap;
            			i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_I_config[ this.level_I_config.length-1].
            				global_haps_string[i].length;j++) {
            			tmp=tmp+ this.level_I_config[ this.level_I_config.length-1].
                				global_haps_string[i][j];
            		}
            		linked_haps.add(tmp);
            	}
            }else {
            	ArrayList<String >  String_X= new ArrayList<String >();
            	ArrayList<String >  String_Y= new ArrayList<String >();
            	ArrayList<String >  String_Z= new ArrayList<String >();
            	int overlap= 1+ this.regions_level_I[2*region_index ][1]- 
            			this.regions_level_II[2*region_index ][0];
//            	
            	
            	for (int i =0;i< this.level_I_config[2*region_index].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_I_config[2*region_index].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_I_config[2*region_index].global_haps_string[i][j];
            		}
            		String_X.add(tmp);
            	}
            	for (int i =0;i< this.level_I_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_I_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_I_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_Y.add(tmp);
            	}
            	
            	for (int i =0;i< this.level_II_config[2*region_index].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_II_config[2*region_index].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_II_config[2*region_index].global_haps_string[i][j];
            		}
            		String_Z.add(tmp);
            	}
//            	
            	for (int i =0;i< String_X.size();i++) {
            		for (int j =0;j< String_Y.size();j++) {
            			if (flag_index==-1) {
            				linked_haps.add(String_X.get(i)+String_Y.get(j) );
    						haps_set.add(String_X.get(i)+String_Y.get(j));
            			}else {
            				if (!haps_set.contains(String_X.get(i)+String_Y.get(j)) ) {
            					haps_set.add(String_X.get(i)+String_Y.get(j));
            					linked_haps.add(String_X.get(i)+String_Y.get(j) );
            				}
	            			for (int k =0;k< String_Z.size();k++) {
	            				String[] comb_str = strmatch( String_X.get(i),
	            						String_Y.get(j), String_Z.get(k),overlap,
	            						flag_index);
	            				for (int h=0;h<comb_str.length;h++ ) {
	            					if (!haps_set.contains(comb_str[h])) {
	            						linked_haps.add(comb_str[h]);
	            						haps_set.add(comb_str[h]);
	            					}
	            				}
	            			}
            			}
            		}
            	}
            }
            
//            for (int i=0;i<linked_haps.size();i++ ) {
//            	System.out.println(linked_haps.get(i).length()); 
//            }
            
            
            
//          
            boolean Using1asDiagonal=true;
            double inpool_site_total_freq=0.0;
            double inpool_site_total_count=0.0;
            for (int l = 0; l < num_site_regional; l++) {
            	for (int i=0; i<inpool_site_freqs[l].length; i++ ) {
            		inpool_site_total_freq+= inpool_site_freqs[l][i];
            		inpool_site_total_count+=1.0;
            	}
            }
            if ((inpool_site_total_freq / inpool_site_total_count) <0.5) {
            	Using1asDiagonal=false;
            }
            System.out.println((inpool_site_total_freq) /inpool_site_total_count);
            
            
            
          //Chen: Construct a diagonal matrix
            
            for (int i=0;i< num_site_regional;i++) {
	            	if (Using1asDiagonal) {
		            	String diagonal_hap =  diagonal_str(i, num_site_regional);
		            	if (!haps_set.contains(diagonal_hap )) {
		            		linked_haps.add(diagonal_hap);
		            	}
	            	}
	            	if (!Using1asDiagonal) {
		            	String diagonal_hap =  diagonal_0_str(i, num_site_regional);
		            	if (!haps_set.contains(diagonal_hap )) {
		            		linked_haps.add(diagonal_hap);
		            	}
	            	}
            }
//            Collections.shuffle(linked_haps);
            
            int haps_2n = linked_haps.size();
            String[][] global_haps_string = new String[haps_2n][num_site_regional];
            String[] hap_IDs = new String[haps_2n];
            for (int h = 0; h < haps_2n; h++) {
                String curr_ID = "";
                String vc_str = linked_haps.get(h);
                String[] vc_arr  = vc_str.split("");
                // the length of vc_str may not reach num_site_regional; so put zeros in.
                for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                    global_haps_string[h][locus] = "0";
                    curr_ID += "0";
                }
                for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
                    global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
                    curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
                }
                hap_IDs[h] = curr_ID;
            }
//            for (int i=0;i<global_haps_string.length;i++ ) {
//            	for (int j=0;j <global_haps_string[i].length;j++ ) {
//            		System.out.print(global_haps_string[i][j]); 
//            	}
//            	System.out.println();
//            }
            //Chen: Construct a diagonal matrix
            
            
            
            
            
            double  diagonal_freq = 0.0002;
            double[] global_haps_freq = new double[haps_2n];
            Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
            
            for (int i=0; i<haps_2n;i++ ) {
            	if ((Using1asDiagonal  ) && (num_of_alternate(global_haps_string[i])==1  )) {
            			global_haps_freq[i]= diagonal_freq;
            	}else if  ((!Using1asDiagonal  ) && (num_of_alternate(global_haps_string[i])
            			==( global_haps_string[i].length-1)  )) {
            		global_haps_freq[i]= diagonal_freq;
            	}
            }

            
            return new HapConfig(
                global_haps_string,
                global_haps_freq,
                null,
                inpool_site_freqs,
                locusInfo,
                this.num_pools,
                hap_IDs,
                pool_IDs,
                this.dp.est_ind_pool,
                Using1asDiagonal);
            
        }
    
    public HapConfig generate_hapconfig_level_IV_fail(
            int[][] the_region, int region_index, String[] pool_IDs, int flag_index) throws IOException {

            int region_start = the_region[region_index][0];
            int region_end = the_region[region_index][1];
            int num_site_regional = region_end-region_start + 1;
            double[][] inpool_site_freqs = new double[num_site_regional][];
            LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];

            for (int l = 0; l < num_site_regional; l++) {
                locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
                inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
            }
            HashSet<String > haps_set =new HashSet<String >();
            ArrayList<String >  linked_haps= new ArrayList<String >();
            if ((this.num_regions_level_I%2==0) && (region_index==this.num_regions_level_IV-1)){
            	for (int i =0;i< this.level_I_config[ this.level_I_config.length-1].num_global_hap;
            			i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_I_config[ this.level_I_config.length-1].
            				global_haps_string[i].length;j++) {
            			tmp=tmp+ this.level_I_config[ this.level_I_config.length-1].
                				global_haps_string[i][j];
            		}
            		linked_haps.add(tmp);
            	}
            }else {
            	ArrayList<String >  String_X= new ArrayList<String >();
            	ArrayList<String >  String_Y= new ArrayList<String >();
            	ArrayList<String >  String_Z= new ArrayList<String >();
            	int overlap= 1+ this.regions_level_I[2*region_index+1 ][1]- 
            			this.regions_level_II[2*region_index +1][0];
            	for (int i =0;i< this.level_I_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_I_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_I_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_X.add(tmp);
            	}
            	for (int i =0;i< this.level_I_config[2*region_index+2].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_I_config[2*region_index+2].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_I_config[2*region_index+2].global_haps_string[i][j];
            		}
            		String_Y.add(tmp);
            	}
            	
            	for (int i =0;i< this.level_II_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_II_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_II_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_Z.add(tmp);
            	}
            	
            	for (int i =0;i< String_X.size();i++) {
            		for (int j =0;j< String_Y.size();j++) {
            			if (flag_index==-1) {
            				linked_haps.add(String_X.get(i)+String_Y.get(j) );
    						haps_set.add(String_X.get(i)+String_Y.get(j));
            			}else {
            				
            				if (!haps_set.contains(String_X.get(i)+String_Y.get(j)) ) {
            					haps_set.add(String_X.get(i)+String_Y.get(j));
            					linked_haps.add(String_X.get(i)+String_Y.get(j) );
            				}
            				for (int k =0;k< String_Z.size();k++) {
            				
	            				String[] comb_str = strmatch( String_X.get(i),
	            						String_Y.get(j), String_Z.get(k),overlap,
	            						flag_index);
	            				for (int h=0;h<comb_str.length;h++ ) {
	            					if (!haps_set.contains(comb_str[h])) {
	            						linked_haps.add(comb_str[h]);
	            						haps_set.add(comb_str[h]);
	            					}
	            				}
                			}
            			}
            		}
            	}
//            	
            	
            }
//            for (int i=0;i<linked_haps.size();i++ ) {
//            	System.out.println(linked_haps.get(i).length()); 
//            }
            
//            
          //Chen: Construct a diagonal matrix
            
            
            boolean Using1asDiagonal=true;
            double inpool_site_total_freq=0.0;
            double inpool_site_total_count=0.0;
            for (int l = 0; l < num_site_regional; l++) {
            	for (int i=0; i<inpool_site_freqs[l].length; i++ ) {
            		inpool_site_total_freq+= inpool_site_freqs[l][i];
            		inpool_site_total_count+=1.0;
            	}
            }
            if ((inpool_site_total_freq / inpool_site_total_count) <0.5) {
            	Using1asDiagonal=false;
            }
            System.out.println((inpool_site_total_freq) /inpool_site_total_count);
            
          //Chen: Construct a diagonal matrix
            
            for (int i=0;i< num_site_regional;i++) {
            	if (Using1asDiagonal) {
	            	String diagonal_hap =  diagonal_str(i, num_site_regional);
	            	if (!haps_set.contains(diagonal_hap )) {
	            		linked_haps.add(diagonal_hap);
	            	}
            	}
            	if (!Using1asDiagonal) {
	            	String diagonal_hap =  diagonal_0_str(i, num_site_regional);
	            	if (!haps_set.contains(diagonal_hap )) {
	            		linked_haps.add(diagonal_hap);
	            	}
            	}
            }
            
            
            
            int haps_2n = linked_haps.size();
            String[][] global_haps_string = new String[haps_2n][num_site_regional];
            String[] hap_IDs = new String[haps_2n];
            for (int h = 0; h < haps_2n; h++) {
                String curr_ID = "";
                String vc_str = linked_haps.get(h);
                String[] vc_arr  = vc_str.split("");
                // the length of vc_str may not reach num_site_regional; so put zeros in.
                for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                    global_haps_string[h][locus] = "0";
                    curr_ID += "0";
                }
                for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
                    global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
                    curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
                }
                hap_IDs[h] = curr_ID;
            }
//            for (int i=0;i<global_haps_string.length;i++ ) {
//            	for (int j=0;j <global_haps_string[i].length;j++ ) {
//            		System.out.print(global_haps_string[i][j]); 
//            	}
//            	System.out.println();
//            }
            //Chen: Construct a diagonal matrix
            
            double[] global_haps_freq = new double[haps_2n];
            Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
            return new HapConfig(
                global_haps_string,
                global_haps_freq,
                null,
                inpool_site_freqs,
                locusInfo,
                this.num_pools,
                hap_IDs,
                pool_IDs,
                this.dp.est_ind_pool,
                Using1asDiagonal);
        }
    
    /**
     *Chen: Construct hap config for level IV
     * 
     */
    
    public HapConfig generate_hapconfig_level_IV(
            int[][] the_region, int region_index, String[] pool_IDs) throws IOException {

            int region_start = the_region[region_index][0];
            int region_end = the_region[region_index][1];
            int num_site_regional = region_end-region_start + 1;
            double[][] inpool_site_freqs = new double[num_site_regional][];
            LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];

            for (int l = 0; l < num_site_regional; l++) {
                locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
                inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
            }
            HashSet<String > haps_set =new HashSet<String >();
            ArrayList<String >  linked_haps= new ArrayList<String >();
            if ((this.num_regions_level_I%2==0) && (region_index==this.num_regions_level_IV-1)){
            	for (int i =0;i< this.level_I_config[ this.level_I_config.length-1].num_global_hap;
            			i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_I_config[ this.level_I_config.length-1].
            				global_haps_string[i].length;j++) {
            			tmp=tmp+ this.level_I_config[ this.level_I_config.length-1].
                				global_haps_string[i][j];
            		}
            		linked_haps.add(tmp);
            	}
            }else {
            	ArrayList<String >  String_X= new ArrayList<String >();
            	ArrayList<String >  String_Y= new ArrayList<String >();
            	ArrayList<String >  String_Z= new ArrayList<String >();
            	int overlap= 1+ this.regions_level_I[2*region_index+1 ][1]- 
            			this.regions_level_II[2*region_index +1][0];
            	for (int i =0;i< this.level_I_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_I_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_I_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_X.add(tmp);
            	}
            	for (int i =0;i< this.level_I_config[2*region_index+2].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_I_config[2*region_index+2].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_I_config[2*region_index+2].global_haps_string[i][j];
            		}
            		String_Y.add(tmp);
            	}
            	
            	for (int i =0;i< this.level_II_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_II_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_II_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_Z.add(tmp);
            	}
            	
            	for (int i =0;i< String_X.size();i++) {
            		for (int j =0;j< String_Y.size();j++) {
            			for (int k =0;k< String_Z.size();k++) {
            				String[] comb_str = strmatch( String_X.get(i),
            						String_Y.get(j), String_Z.get(k),overlap,
            						this.level_III_IV_region_mismatch_tolerance);
            				if   ((comb_str.length>0) &&
            						(!haps_set.contains(String_X.get(i)+  String_Y.get(j)))) {
            					linked_haps.add(String_X.get(i)+  String_Y.get(j));
        						haps_set.add(String_X.get(i)+  String_Y.get(j));
        						break;
            				}
            				for (int h=0;h<comb_str.length;h++ ) {
            					if (!haps_set.contains(comb_str[h])) {
            						linked_haps.add(comb_str[h]);
            						haps_set.add(comb_str[h]);
            					}
            				}
            			}
            		}
            	}
//            	
            	
            }
//            for (int i=0;i<linked_haps.size();i++ ) {
//            	System.out.println(linked_haps.get(i).length()); 
//            }
            
//            
          //Chen: Construct a diagonal matrix
            
            
            boolean Using1asDiagonal=true;
            double inpool_site_total_freq=0.0;
            double inpool_site_total_count=0.0;
            for (int l = 0; l < num_site_regional; l++) {
            	for (int i=0; i<inpool_site_freqs[l].length; i++ ) {
            		inpool_site_total_freq+= inpool_site_freqs[l][i];
            		inpool_site_total_count+=1.0;
            	}
            }
            if ((inpool_site_total_freq / inpool_site_total_count) <0.5) {
            	Using1asDiagonal=false;
            }
            System.out.println((inpool_site_total_freq) /inpool_site_total_count);
            
          //Chen: Construct a diagonal matrix
            
            for (int i=0;i< num_site_regional;i++) {
//            	if (Using1asDiagonal) {
//	            	String diagonal_hap =  diagonal_str(i, num_site_regional);
//	            	if (!haps_set.contains(diagonal_hap )) {
//	            		linked_haps.add(diagonal_hap);
//	            	}
//            	}
//            	if (!Using1asDiagonal) {
//	            	String diagonal_hap =  diagonal_0_str(i, num_site_regional);
//	            	if (!haps_set.contains(diagonal_hap )) {
//	            		linked_haps.add(diagonal_hap);
//	            	}
//            	}
            }
            
            
            
            int haps_2n = linked_haps.size();
            String[][] global_haps_string = new String[haps_2n][num_site_regional];
            String[] hap_IDs = new String[haps_2n];
            for (int h = 0; h < haps_2n; h++) {
                String curr_ID = "";
                String vc_str = linked_haps.get(h);
                String[] vc_arr  = vc_str.split("");
                // the length of vc_str may not reach num_site_regional; so put zeros in.
                for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                    global_haps_string[h][locus] = "0";
                    curr_ID += "0";
                }
                for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
                    global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
                    curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
                }
                hap_IDs[h] = curr_ID;
            }
//            for (int i=0;i<global_haps_string.length;i++ ) {
//            	for (int j=0;j <global_haps_string[i].length;j++ ) {
//            		System.out.print(global_haps_string[i][j]); 
//            	}
//            	System.out.println();
//            }
            //Chen: Construct a diagonal matrix
            
            double[] global_haps_freq = new double[haps_2n];
            Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
            return new HapConfig(
                global_haps_string,
                global_haps_freq,
                null,
                inpool_site_freqs,
                locusInfo,
                this.num_pools,
                hap_IDs,
                pool_IDs,
                this.dp.est_ind_pool,
                Using1asDiagonal);
        }
    
    
    public HapConfig generate_hapconfig_level_V_fail(
            int[][] the_region, int region_index, String[] pool_IDs, int flag_index) throws IOException {

            int region_start = the_region[region_index][0];
            int region_end = the_region[region_index][1];
            int num_site_regional = region_end-region_start + 1;
            double[][] inpool_site_freqs = new double[num_site_regional][];
            LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];

            for (int l = 0; l < num_site_regional; l++) {
                locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
                inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
            }
            HashSet<String > haps_set =new HashSet<String >();
            ArrayList<String >  linked_haps= new ArrayList<String >();
            if ((this.num_regions_level_III%2==1) && (region_index==this.num_regions_level_V-1)){
            	for (int i =0;i< this.level_III_config[ this.level_III_config.length-1].num_global_hap;
            			i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_III_config[ this.level_III_config.length-1].
            				global_haps_string[i].length;j++) {
            			tmp=tmp+ this.level_III_config[ this.level_III_config.length-1].
                				global_haps_string[i][j];
            		}
            		linked_haps.add(tmp);
            	}
            }else {
            	ArrayList<String >  String_X= new ArrayList<String >();
            	ArrayList<String >  String_Y= new ArrayList<String >();
            	ArrayList<String >  String_Z= new ArrayList<String >();
            	int overlap= 1+ this.regions_level_III[2*region_index ][1]- 
            			this.regions_level_IV[2*region_index ][0];
//            	
            	for (int i =0;i< this.level_III_config[2*region_index].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_III_config[2*region_index].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_III_config[2*region_index].global_haps_string[i][j];
            		}
            		String_X.add(tmp);
            	}
            	for (int i =0;i< this.level_III_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_III_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_III_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_Y.add(tmp);
            	}
            	
            	for (int i =0;i< this.level_IV_config[2*region_index].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_IV_config[2*region_index].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_IV_config[2*region_index].global_haps_string[i][j];
            		}
            		String_Z.add(tmp);
            	}
            	
            	
//            	
            	for (int i =0;i< String_X.size();i++) {
            		for (int j =0;j< String_Y.size();j++) {
            			if (flag_index==-1) {
            				linked_haps.add(String_X.get(i)+String_Y.get(j) );
        					haps_set.add(String_X.get(i)+String_Y.get(j));
            			}else  if (flag_index!=-2){
            				if (!haps_set.contains(String_X.get(i)+String_Y.get(j)) ) {
            					haps_set.add(String_X.get(i)+String_Y.get(j));
            					linked_haps.add(String_X.get(i)+String_Y.get(j) );
            				}
            				for (int k =0;k< String_Z.size();k++) {
                				String[] comb_str = strmatch( String_X.get(i),
                						String_Y.get(j), String_Z.get(k),overlap,
                						flag_index);
                				for (int h=0;h<comb_str.length;h++ ) {
                					if (!haps_set.contains(comb_str[h])) {
                						linked_haps.add(comb_str[h]);
                						haps_set.add(comb_str[h]);
                					}
                				}
                			}
            			}else  if (flag_index==-2){
            				for (int k =0;k< String_Z.size();k++) {
	            				String[] comb_str = strmatch( String_X.get(i),
	            						String_Y.get(j), String_Z.get(k),overlap,
	            						this.level_V_VI_region_mismatch_tolerance);
	            				if (comb_str.length!=0) {
	            					String tmp_hap= String_X.get(i)+String_Y.get(j);
	            					if (!haps_set.contains(tmp_hap) ) {
	                					haps_set.add(tmp_hap);
	                					linked_haps.add(tmp_hap);
	                				}
	            					
		            				tmp_hap = tmp_hap.substring(0,this.regions_level_IV[2*region_index ][0]
		            						- this.regions_level_III[2*region_index ][0])
		            				+String_Z.get(k)
									+tmp_hap.substring(this.regions_level_IV[2*region_index ][0]
											-this.regions_level_III[2*region_index ][0]
											+ String_Z.get(k).length());
		            				if (!haps_set.contains(tmp_hap) ) {
	                					haps_set.add(tmp_hap);
	                					linked_haps.add(tmp_hap);
	                				}
	            				}
            				}
            			}
            		}
            	}
            }
 
//            
          //Chen: Construct a diagonal matrix
            
            boolean Using1asDiagonal=true;
            double inpool_site_total_freq=0.0;
            double inpool_site_total_count=0.0;
            for (int l = 0; l < num_site_regional; l++) {
            	for (int i=0; i<inpool_site_freqs[l].length; i++ ) {
            		inpool_site_total_freq+= inpool_site_freqs[l][i];
            		inpool_site_total_count+=1.0;
            	}
            }
            if ((inpool_site_total_freq / inpool_site_total_count) <0.5) {
            	Using1asDiagonal=false;
            }
            System.out.println((inpool_site_total_freq) /inpool_site_total_count);
            
          //Chen: Construct a diagonal matrix
            
            for (int i=0;i< num_site_regional;i++) {
            	if (Using1asDiagonal) {
	            	String diagonal_hap =  diagonal_str(i, num_site_regional);
	            	if (!haps_set.contains(diagonal_hap )) {
	            		linked_haps.add(diagonal_hap);
	            	}
            	}
            	if (!Using1asDiagonal) {
	            	String diagonal_hap =  diagonal_0_str(i, num_site_regional);
	            	if (!haps_set.contains(diagonal_hap )) {
	            		linked_haps.add(diagonal_hap);
	            	}
            	}
            }
            
            
            
            int haps_2n = linked_haps.size();
            String[][] global_haps_string = new String[haps_2n][num_site_regional];
            String[] hap_IDs = new String[haps_2n];
            for (int h = 0; h < haps_2n; h++) {
                String curr_ID = "";
                String vc_str = linked_haps.get(h);
                String[] vc_arr  = vc_str.split("");
                // the length of vc_str may not reach num_site_regional; so put zeros in.
                for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                    global_haps_string[h][locus] = "0";
                    curr_ID += "0";
                }
                for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
                    global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
                    curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
                }
                hap_IDs[h] = curr_ID;
            }
//            for (int i=0;i<global_haps_string.length;i++ ) {
//            	for (int j=0;j <global_haps_string[i].length;j++ ) {
//            		System.out.print(global_haps_string[i][j]); 
//            	}
//            	System.out.println();
//            }
            //Chen: Construct a diagonal matrix
            
            double[] global_haps_freq = new double[haps_2n];
            Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
            return new HapConfig(
                global_haps_string,
                global_haps_freq,
                null,
                inpool_site_freqs,
                locusInfo,
                this.num_pools,
                hap_IDs,
                pool_IDs,
                this.dp.est_ind_pool,
                Using1asDiagonal);
        }
    
    
    
    
    
    public HapConfig generate_hapconfig_level_V(
            int[][] the_region, int region_index, String[] pool_IDs) throws IOException {

            int region_start = the_region[region_index][0];
            int region_end = the_region[region_index][1];
            int num_site_regional = region_end-region_start + 1;
            double[][] inpool_site_freqs = new double[num_site_regional][];
            LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];

            for (int l = 0; l < num_site_regional; l++) {
                locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
                inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
            }
            HashSet<String > haps_set =new HashSet<String >();
            ArrayList<String >  linked_haps= new ArrayList<String >();
            if ((this.num_regions_level_III%2==1) && (region_index==this.num_regions_level_V-1)){
            	for (int i =0;i< this.level_III_config[ this.level_III_config.length-1].num_global_hap;
            			i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_III_config[ this.level_III_config.length-1].
            				global_haps_string[i].length;j++) {
            			tmp=tmp+ this.level_III_config[ this.level_III_config.length-1].
                				global_haps_string[i][j];
            		}
            		linked_haps.add(tmp);
            	}
            }else {
            	ArrayList<String >  String_X= new ArrayList<String >();
            	ArrayList<String >  String_Y= new ArrayList<String >();
            	ArrayList<String >  String_Z= new ArrayList<String >();
            	int overlap= 1+ this.regions_level_III[2*region_index ][1]- 
            			this.regions_level_IV[2*region_index ][0];
//            	
            	for (int i =0;i< this.level_III_config[2*region_index].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_III_config[2*region_index].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_III_config[2*region_index].global_haps_string[i][j];
            		}
            		String_X.add(tmp);
            	}
            	for (int i =0;i< this.level_III_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_III_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_III_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_Y.add(tmp);
            	}
            	
            	for (int i =0;i< this.level_IV_config[2*region_index].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_IV_config[2*region_index].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_IV_config[2*region_index].global_haps_string[i][j];
            		}
            		String_Z.add(tmp);
            	}
//            	
            	for (int i =0;i< String_X.size();i++) {
            		for (int j =0;j< String_Y.size();j++) {
            			for (int k =0;k< String_Z.size();k++) {
            				String[] comb_str = strmatch( String_X.get(i),
            						String_Y.get(j), String_Z.get(k),overlap,
            						this.level_V_VI_region_mismatch_tolerance);
            				if   ((comb_str.length>0) &&
            						(!haps_set.contains(String_X.get(i)+  String_Y.get(j)))) {
            					linked_haps.add(String_X.get(i)+  String_Y.get(j));
        						haps_set.add(String_X.get(i)+  String_Y.get(j));
        						break;
            				}
            				for (int h=0;h<comb_str.length;h++ ) {
            					if (!haps_set.contains(comb_str[h])) {
            						linked_haps.add(comb_str[h]);
            						haps_set.add(comb_str[h]);
            					}
            				}
            			}
            		}
            	}
            }
 
//            
          //Chen: Construct a diagonal matrix
            
            boolean Using1asDiagonal=true;
            double inpool_site_total_freq=0.0;
            double inpool_site_total_count=0.0;
            for (int l = 0; l < num_site_regional; l++) {
            	for (int i=0; i<inpool_site_freqs[l].length; i++ ) {
            		inpool_site_total_freq+= inpool_site_freqs[l][i];
            		inpool_site_total_count+=1.0;
            	}
            }
            if ((inpool_site_total_freq / inpool_site_total_count) <0.5) {
            	Using1asDiagonal=false;
            }
            System.out.println((inpool_site_total_freq) /inpool_site_total_count);
            
          //Chen: Construct a diagonal matrix
            
            for (int i=0;i< num_site_regional;i++) {
//            	if (Using1asDiagonal) {
//	            	String diagonal_hap =  diagonal_str(i, num_site_regional);
//	            	if (!haps_set.contains(diagonal_hap )) {
//	            		linked_haps.add(diagonal_hap);
//	            	}
//            	}
//            	if (!Using1asDiagonal) {
//	            	String diagonal_hap =  diagonal_0_str(i, num_site_regional);
//	            	if (!haps_set.contains(diagonal_hap )) {
//	            		linked_haps.add(diagonal_hap);
//	            	}
//            	}
            }
            
            
            
            int haps_2n = linked_haps.size();
            String[][] global_haps_string = new String[haps_2n][num_site_regional];
            String[] hap_IDs = new String[haps_2n];
            for (int h = 0; h < haps_2n; h++) {
                String curr_ID = "";
                String vc_str = linked_haps.get(h);
                String[] vc_arr  = vc_str.split("");
                // the length of vc_str may not reach num_site_regional; so put zeros in.
                for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                    global_haps_string[h][locus] = "0";
                    curr_ID += "0";
                }
                for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
                    global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
                    curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
                }
                hap_IDs[h] = curr_ID;
            }
//            for (int i=0;i<global_haps_string.length;i++ ) {
//            	for (int j=0;j <global_haps_string[i].length;j++ ) {
//            		System.out.print(global_haps_string[i][j]); 
//            	}
//            	System.out.println();
//            }
            //Chen: Construct a diagonal matrix
            
            double[] global_haps_freq = new double[haps_2n];
            Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
            return new HapConfig(
                global_haps_string,
                global_haps_freq,
                null,
                inpool_site_freqs,
                locusInfo,
                this.num_pools,
                hap_IDs,
                pool_IDs,
                this.dp.est_ind_pool,
                Using1asDiagonal);
        }
    
    
    
    public HapConfig generate_hapconfig_level_VII(
            int[][] the_region, int region_index, String[] pool_IDs) throws IOException {

            int region_start = the_region[region_index][0];
            int region_end = the_region[region_index][1];
            int num_site_regional = region_end-region_start + 1;
            double[][] inpool_site_freqs = new double[num_site_regional][];
            LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];

            for (int l = 0; l < num_site_regional; l++) {
                locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
                inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
            }
            HashSet<String > haps_set =new HashSet<String >();
            ArrayList<String >  linked_haps= new ArrayList<String >();
            if ((this.num_regions_level_V%2==1) && (region_index==this.num_regions_level_VII-1)){
            	for (int i =0;i< this.level_V_config[ this.level_V_config.length-1].num_global_hap;
            			i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_V_config[ this.level_V_config.length-1].
            				global_haps_string[i].length;j++) {
            			tmp=tmp+ this.level_V_config[ this.level_V_config.length-1].
                				global_haps_string[i][j];
            		}
            		linked_haps.add(tmp);
            	}
            }else {
            	ArrayList<String >  String_X= new ArrayList<String >();
            	ArrayList<String >  String_Y= new ArrayList<String >();
            	ArrayList<String >  String_Z= new ArrayList<String >();
            	int overlap= 1+ this.regions_level_V[2*region_index ][1]- 
            			this.regions_level_VI[2*region_index ][0];
//            	
            	for (int i =0;i< this.level_V_config[2*region_index].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_V_config[2*region_index].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_V_config[2*region_index].global_haps_string[i][j];
            		}
            		String_X.add(tmp);
            	}
            	for (int i =0;i< this.level_V_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_V_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_V_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_Y.add(tmp);
            	}
            	
            	for (int i =0;i< this.level_VI_config[2*region_index].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_VI_config[2*region_index].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_VI_config[2*region_index].global_haps_string[i][j];
            		}
            		String_Z.add(tmp);
            	}
//            	
            	for (int i =0;i< String_X.size();i++) {
            		for (int j =0;j< String_Y.size();j++) {
            			for (int k =0;k< String_Z.size();k++) {
            				String[] comb_str = strmatch( String_X.get(i),
            						String_Y.get(j), String_Z.get(k),overlap,
            						this.level_VII_VIII_region_mismatch_tolerance);
            				if   ((comb_str.length>0) &&
            						(!haps_set.contains(String_X.get(i)+  String_Y.get(j)))) {
            					linked_haps.add(String_X.get(i)+  String_Y.get(j));
        						haps_set.add(String_X.get(i)+  String_Y.get(j));
        						break;
            				}
            				for (int h=0;h<comb_str.length;h++ ) {
            					if (!haps_set.contains(comb_str[h])) {
            						linked_haps.add(comb_str[h]);
            						haps_set.add(comb_str[h]);
            					}
            				}
            			}
            		}
            	}
            }
 
//            
          //Chen: Construct a diagonal matrix
            
            boolean Using1asDiagonal=true;
            double inpool_site_total_freq=0.0;
            double inpool_site_total_count=0.0;
            for (int l = 0; l < num_site_regional; l++) {
            	for (int i=0; i<inpool_site_freqs[l].length; i++ ) {
            		inpool_site_total_freq+= inpool_site_freqs[l][i];
            		inpool_site_total_count+=1.0;
            	}
            }
            if ((inpool_site_total_freq / inpool_site_total_count) <0.5) {
            	Using1asDiagonal=false;
            }
            System.out.println((inpool_site_total_freq) /inpool_site_total_count);
            
          //Chen: Construct a diagonal matrix
            
            for (int i=0;i< num_site_regional;i++) {
//            	if (Using1asDiagonal) {
//	            	String diagonal_hap =  diagonal_str(i, num_site_regional);
//	            	if (!haps_set.contains(diagonal_hap )) {
//	            		linked_haps.add(diagonal_hap);
//	            	}
//            	}
//            	if (!Using1asDiagonal) {
//	            	String diagonal_hap =  diagonal_0_str(i, num_site_regional);
//	            	if (!haps_set.contains(diagonal_hap )) {
//	            		linked_haps.add(diagonal_hap);
//	            	}
//            	}
            }
            
            
            
            int haps_2n = linked_haps.size();
            String[][] global_haps_string = new String[haps_2n][num_site_regional];
            String[] hap_IDs = new String[haps_2n];
            for (int h = 0; h < haps_2n; h++) {
                String curr_ID = "";
                String vc_str = linked_haps.get(h);
                String[] vc_arr  = vc_str.split("");
                // the length of vc_str may not reach num_site_regional; so put zeros in.
                for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                    global_haps_string[h][locus] = "0";
                    curr_ID += "0";
                }
                for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
                    global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
                    curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
                }
                hap_IDs[h] = curr_ID;
            }
//            for (int i=0;i<global_haps_string.length;i++ ) {
//            	for (int j=0;j <global_haps_string[i].length;j++ ) {
//            		System.out.print(global_haps_string[i][j]); 
//            	}
//            	System.out.println();
//            }
            //Chen: Construct a diagonal matrix
            
            double[] global_haps_freq = new double[haps_2n];
            Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
            return new HapConfig(
                global_haps_string,
                global_haps_freq,
                null,
                inpool_site_freqs,
                locusInfo,
                this.num_pools,
                hap_IDs,
                pool_IDs,
                this.dp.est_ind_pool,
                Using1asDiagonal);
        }
    
    
    
    public HapConfig generate_hapconfig_level_VII_fail(
    		int[][] the_region, int region_index, String[] pool_IDs, int flag_index) throws IOException {

            int region_start = the_region[region_index][0];
            int region_end = the_region[region_index][1];
            int num_site_regional = region_end-region_start + 1;
            double[][] inpool_site_freqs = new double[num_site_regional][];
            LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];

            for (int l = 0; l < num_site_regional; l++) {
                locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
                inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
            }
            HashSet<String > haps_set =new HashSet<String >();
            ArrayList<String >  linked_haps= new ArrayList<String >();
            if ((this.num_regions_level_V%2==1) && (region_index==this.num_regions_level_VII-1)){
            	for (int i =0;i< this.level_V_config[ this.level_V_config.length-1].num_global_hap;
            			i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_V_config[ this.level_V_config.length-1].
            				global_haps_string[i].length;j++) {
            			tmp=tmp+ this.level_V_config[ this.level_V_config.length-1].
                				global_haps_string[i][j];
            		}
            		linked_haps.add(tmp);
            	}
            }else {
            	ArrayList<String >  String_X= new ArrayList<String >();
            	ArrayList<String >  String_Y= new ArrayList<String >();
            	ArrayList<String >  String_Z= new ArrayList<String >();
            	int overlap= 1+ this.regions_level_V[2*region_index ][1]- 
            			this.regions_level_VI[2*region_index ][0];
//            	
            	for (int i =0;i< this.level_V_config[2*region_index].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_V_config[2*region_index].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_V_config[2*region_index].global_haps_string[i][j];
            		}
            		String_X.add(tmp);
            	}
            	for (int i =0;i< this.level_V_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_V_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_V_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_Y.add(tmp);
            	}
            	
            	for (int i =0;i< this.level_VI_config[2*region_index].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_VI_config[2*region_index].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_VI_config[2*region_index].global_haps_string[i][j];
            		}
            		String_Z.add(tmp);
            	}
//            	
            	for (int i =0;i< String_X.size();i++) {
            		for (int j =0;j< String_Y.size();j++) {
            			if (flag_index==-1) {
            				linked_haps.add(String_X.get(i)+String_Y.get(j) );
        					haps_set.add(String_X.get(i)+String_Y.get(j));
            			}else  if (flag_index!=-2){
            				if (!haps_set.contains(String_X.get(i)+String_Y.get(j)) ) {
            					haps_set.add(String_X.get(i)+String_Y.get(j));
            					linked_haps.add(String_X.get(i)+String_Y.get(j) );
            				}
            				for (int k =0;k< String_Z.size();k++) {
                				String[] comb_str = strmatch( String_X.get(i),
                						String_Y.get(j), String_Z.get(k),overlap,
                						flag_index);
                				for (int h=0;h<comb_str.length;h++ ) {
                					if (!haps_set.contains(comb_str[h])) {
                						linked_haps.add(comb_str[h]);
                						haps_set.add(comb_str[h]);
                					}
                				}
                			}
            			}else  if (flag_index==-2){
            				for (int k =0;k< String_Z.size();k++) {
	            				String[] comb_str = strmatch( String_X.get(i),
	            						String_Y.get(j), String_Z.get(k),overlap,
	            						this.level_VII_VIII_region_mismatch_tolerance);
	            				if (comb_str.length!=0) {
	            					String tmp_hap= String_X.get(i)+String_Y.get(j);
	            					if (!haps_set.contains(tmp_hap) ) {
	                					haps_set.add(tmp_hap);
	                					linked_haps.add(tmp_hap);
	                					
	                				}
	            					
		            				tmp_hap = tmp_hap.substring(0,this.regions_level_VI[2*region_index][0]-
		            						this.regions_level_V[2*region_index][0])
		            				+String_Z.get(k)
									+tmp_hap.substring(this.regions_level_VI[2*region_index][0]-
											this.regions_level_V[2*region_index][0]
											+ String_Z.get(k).length());
		            				if (!haps_set.contains(tmp_hap) ) {
	                					haps_set.add(tmp_hap);
	                					linked_haps.add(tmp_hap);
	                				}
	            				}
            				}
            			}
            		}
            	}
            }
 
//            
          //Chen: Construct a diagonal matrix
            
            boolean Using1asDiagonal=true;
            double inpool_site_total_freq=0.0;
            double inpool_site_total_count=0.0;
            for (int l = 0; l < num_site_regional; l++) {
            	for (int i=0; i<inpool_site_freqs[l].length; i++ ) {
            		inpool_site_total_freq+= inpool_site_freqs[l][i];
            		inpool_site_total_count+=1.0;
            	}
            }
            if ((inpool_site_total_freq / inpool_site_total_count) <0.5) {
            	Using1asDiagonal=false;
            }
            System.out.println((inpool_site_total_freq) /inpool_site_total_count);
            
          //Chen: Construct a diagonal matrix
            
            for (int i=0;i< num_site_regional;i++) {
            	if (Using1asDiagonal) {
	            	String diagonal_hap =  diagonal_str(i, num_site_regional);
	            	if (!haps_set.contains(diagonal_hap )) {
	            		linked_haps.add(diagonal_hap);
	            	}
            	}
            	if (!Using1asDiagonal) {
	            	String diagonal_hap =  diagonal_0_str(i, num_site_regional);
	            	if (!haps_set.contains(diagonal_hap )) {
	            		linked_haps.add(diagonal_hap);
	            	}
            	}
            }
            
            
            
            int haps_2n = linked_haps.size();
            String[][] global_haps_string = new String[haps_2n][num_site_regional];
            String[] hap_IDs = new String[haps_2n];
            for (int h = 0; h < haps_2n; h++) {
                String curr_ID = "";
                String vc_str = linked_haps.get(h);
                String[] vc_arr  = vc_str.split("");
                // the length of vc_str may not reach num_site_regional; so put zeros in.
                for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                    global_haps_string[h][locus] = "0";
                    curr_ID += "0";
                }
                for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
                    global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
                    curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
                }
                hap_IDs[h] = curr_ID;
            }
//            for (int i=0;i<global_haps_string.length;i++ ) {
//            	for (int j=0;j <global_haps_string[i].length;j++ ) {
//            		System.out.print(global_haps_string[i][j]); 
//            	}
//            	System.out.println();
//            }
            //Chen: Construct a diagonal matrix
            
            double[] global_haps_freq = new double[haps_2n];
            Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
            return new HapConfig(
                global_haps_string,
                global_haps_freq,
                null,
                inpool_site_freqs,
                locusInfo,
                this.num_pools,
                hap_IDs,
                pool_IDs,
                this.dp.est_ind_pool,
                Using1asDiagonal);
    }
    
    
    public HapConfig generate_hapconfig_level_VI_fail(
	            int[][] the_region, int region_index, String[] pool_IDs, int flag_index
	            ) throws IOException {
	
	            int region_start = the_region[region_index][0];
	            int region_end = the_region[region_index][1];
	            int num_site_regional = region_end-region_start + 1;
	            double[][] inpool_site_freqs = new double[num_site_regional][];
	            LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];
	
	            for (int l = 0; l < num_site_regional; l++) {
	                locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
	                inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
	            }
	            HashSet<String > haps_set =new HashSet<String >();
	            ArrayList<String >  linked_haps= new ArrayList<String >();
	            if ((this.num_regions_level_III%2==0) && (region_index==this.num_regions_level_VI-1)){
	            	for (int i =0;i< this.level_III_config[ this.level_III_config.length-1].num_global_hap;
	            			i++) {
	            		String tmp ="";
	            		for (int j =0; j< this.level_III_config[ this.level_III_config.length-1].
	            				global_haps_string[i].length;j++) {
	            			tmp=tmp+ this.level_III_config[ this.level_III_config.length-1].
	                				global_haps_string[i][j];
	            		}
	            		linked_haps.add(tmp);
	            	}
	            }else {
	            	ArrayList<String >  String_X= new ArrayList<String >();
	            	ArrayList<String >  String_Y= new ArrayList<String >();
	            	ArrayList<String >  String_Z= new ArrayList<String >();
	            	int overlap= 1+ this.regions_level_III[2*region_index+1 ][1]- 
	            			this.regions_level_IV[2*region_index +1][0];
	//            	
	            	int mismatch_cutoff=0;
	            	for (int i =0;i< this.level_III_config[2*region_index+1].num_global_hap;i++) {
	            		String tmp ="";
	            		for (int j =0; j< this.level_III_config[2*region_index+1].global_haps_string[i].length;
	            				j++) {
	            			tmp=tmp+ this.level_III_config[2*region_index+1].global_haps_string[i][j];
	            		}
	            		String_X.add(tmp);
	            	}
	            	for (int i =0;i< this.level_III_config[2*region_index+2].num_global_hap;i++) {
	            		String tmp ="";
	            		for (int j =0; j< this.level_III_config[2*region_index+2].global_haps_string[i].length;
	            				j++) {
	            			tmp=tmp+ this.level_III_config[2*region_index+2].global_haps_string[i][j];
	            		}
	            		String_Y.add(tmp);
	            	}
	            	
	            	for (int i =0;i< this.level_IV_config[2*region_index+1].num_global_hap;i++) {
	            		String tmp ="";
	            		for (int j =0; j< this.level_IV_config[2*region_index+1].global_haps_string[i].length;
	            				j++) {
	            			tmp=tmp+ this.level_IV_config[2*region_index+1].global_haps_string[i][j];
	            		}
	            		String_Z.add(tmp);
	            	}
	//            	
	            	for (int i =0;i< String_X.size();i++) {
	            		for (int j =0;j< String_Y.size();j++) {
	            			if (flag_index==-1) {
	            				linked_haps.add(String_X.get(i)+String_Y.get(j) );
	        					haps_set.add(String_X.get(i)+String_Y.get(j));
	            			}else  if (flag_index!=-2){
	            				if (!haps_set.contains(String_X.get(i)+String_Y.get(j)) ) {
	            					haps_set.add(String_X.get(i)+String_Y.get(j));
	            					linked_haps.add(String_X.get(i)+String_Y.get(j) );
	            				}
	            				for (int k =0;k< String_Z.size();k++) {
	                				String[] comb_str = strmatch( String_X.get(i),
	                						String_Y.get(j), String_Z.get(k),overlap,
	                						flag_index);
	                				for (int h=0;h<comb_str.length;h++ ) {
	                					if (!haps_set.contains(comb_str[h])) {
	                						linked_haps.add(comb_str[h]);
	                						haps_set.add(comb_str[h]);
	                					}
	                				}
	                			}
	            			}else  if (flag_index==-2){
	            				for (int k =0;k< String_Z.size();k++) {
		            				String[] comb_str = strmatch( String_X.get(i),
		            						String_Y.get(j), String_Z.get(k),overlap,
		            						this.level_V_VI_region_mismatch_tolerance);
		            				if (comb_str.length!=0) {
		            					String tmp_hap= String_X.get(i)+String_Y.get(j);
		            					if (!haps_set.contains(tmp_hap) ) {
		                					haps_set.add(tmp_hap);
		                					linked_haps.add(tmp_hap);
		                				}
		            					
			            				tmp_hap = tmp_hap.substring(0,this.regions_level_IV[2*region_index+1][0]- 
			            						this.regions_level_III[2*region_index+1 ][0])
			            				+String_Z.get(k)
										+tmp_hap.substring(this.regions_level_IV[2*region_index+1][0]-
												this.regions_level_III[2*region_index+1 ][0]
											+ String_Z.get(k).length());
			            				if (!haps_set.contains(tmp_hap) ) {
		                					haps_set.add(tmp_hap);
		                					linked_haps.add(tmp_hap);
		                				}
		            				}
	            				}
	            			}
	            		}
	            	}
	            }
	//            for (int i=0;i<linked_haps.size();i++ ) {
	//            	System.out.println(linked_haps.get(i).length()); 
	//            }
	            
	//            
	          //Chen: Construct a diagonal matrix
	            
	            
	            boolean Using1asDiagonal=true;
	            double inpool_site_total_freq=0.0;
	            double inpool_site_total_count=0.0;
	            for (int l = 0; l < num_site_regional; l++) {
	            	for (int i=0; i<inpool_site_freqs[l].length; i++ ) {
	            		inpool_site_total_freq+= inpool_site_freqs[l][i];
	            		inpool_site_total_count+=1.0;
	            	}
	            }
	            if ((inpool_site_total_freq / inpool_site_total_count) <0.5) {
	            	Using1asDiagonal=false;
	            }
	            System.out.println((inpool_site_total_freq) /inpool_site_total_count);
	            
	          //Chen: Construct a diagonal matrix
	            
	            for (int i=0;i< num_site_regional;i++) {
	            	if (Using1asDiagonal) {
		            	String diagonal_hap =  diagonal_str(i, num_site_regional);
		            	if (!haps_set.contains(diagonal_hap )) {
		            		linked_haps.add(diagonal_hap);
		            	}
	            	}
	            	if (!Using1asDiagonal) {
		            	String diagonal_hap =  diagonal_0_str(i, num_site_regional);
		            	if (!haps_set.contains(diagonal_hap )) {
		            		linked_haps.add(diagonal_hap);
		            	}
	            	}
	            }
	            
	            int haps_2n = linked_haps.size();
	            String[][] global_haps_string = new String[haps_2n][num_site_regional];
	            String[] hap_IDs = new String[haps_2n];
	            for (int h = 0; h < haps_2n; h++) {
	                String curr_ID = "";
	                String vc_str = linked_haps.get(h);
	                String[] vc_arr  = vc_str.split("");
	                // the length of vc_str may not reach num_site_regional; so put zeros in.
	                for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
	                    global_haps_string[h][locus] = "0";
	                    curr_ID += "0";
	                }
	                for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
	                    global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
	                    curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
	                }
	                hap_IDs[h] = curr_ID;
	            }
	//            for (int i=0;i<global_haps_string.length;i++ ) {
	//            	for (int j=0;j <global_haps_string[i].length;j++ ) {
	//            		System.out.print(global_haps_string[i][j]); 
	//            	}
	//            	System.out.println();
	//            }
	            //Chen: Construct a diagonal matrix
	            
	            double[] global_haps_freq = new double[haps_2n];
	            Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
	            return new HapConfig(
	                global_haps_string,
	                global_haps_freq,
	                null,
	                inpool_site_freqs,
	                locusInfo,
	                this.num_pools,
	                hap_IDs,
	                pool_IDs,
	                this.dp.est_ind_pool,
	                Using1asDiagonal);
	        }

	public HapConfig generate_hapconfig_level_VI(
            int[][] the_region, int region_index, String[] pool_IDs) throws IOException {

            int region_start = the_region[region_index][0];
            int region_end = the_region[region_index][1];
            int num_site_regional = region_end-region_start + 1;
            double[][] inpool_site_freqs = new double[num_site_regional][];
            LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];

            for (int l = 0; l < num_site_regional; l++) {
                locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
                inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
            }
            HashSet<String > haps_set =new HashSet<String >();
            ArrayList<String >  linked_haps= new ArrayList<String >();
            if ((this.num_regions_level_III%2==0) && (region_index==this.num_regions_level_VI-1)){
            	for (int i =0;i< this.level_III_config[ this.level_III_config.length-1].num_global_hap;
            			i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_III_config[ this.level_III_config.length-1].
            				global_haps_string[i].length;j++) {
            			tmp=tmp+ this.level_III_config[ this.level_III_config.length-1].
                				global_haps_string[i][j];
            		}
            		linked_haps.add(tmp);
            	}
            }else {
            	ArrayList<String >  String_X= new ArrayList<String >();
            	ArrayList<String >  String_Y= new ArrayList<String >();
            	ArrayList<String >  String_Z= new ArrayList<String >();
            	int overlap= 1+ this.regions_level_III[2*region_index+1 ][1]- 
            			this.regions_level_IV[2*region_index +1][0];
//            	
            	int mismatch_cutoff=0;
            	for (int i =0;i< this.level_III_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_III_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_III_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_X.add(tmp);
            	}
            	for (int i =0;i< this.level_III_config[2*region_index+2].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_III_config[2*region_index+2].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_III_config[2*region_index+2].global_haps_string[i][j];
            		}
            		String_Y.add(tmp);
            	}
            	
            	for (int i =0;i< this.level_IV_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_IV_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_IV_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_Z.add(tmp);
            	}
//            	
            	for (int i =0;i< String_X.size();i++) {
            		for (int j =0;j< String_Y.size();j++) {
            			for (int k =0;k< String_Z.size();k++) {
            				String[] comb_str = strmatch( String_X.get(i),
            						String_Y.get(j), String_Z.get(k),overlap,
            						this.level_V_VI_region_mismatch_tolerance);
            				if   ((comb_str.length>0) &&
            						(!haps_set.contains(String_X.get(i)+  String_Y.get(j)))) {
            					linked_haps.add(String_X.get(i)+  String_Y.get(j));
        						haps_set.add(String_X.get(i)+  String_Y.get(j));
        						break;
            				}
            				for (int h=0;h<comb_str.length;h++ ) {
            					if (!haps_set.contains(comb_str[h])) {
            						linked_haps.add(comb_str[h]);
            						haps_set.add(comb_str[h]);
            					}
            				}
            			}
            		}
            	}
            }
//            for (int i=0;i<linked_haps.size();i++ ) {
//            	System.out.println(linked_haps.get(i).length()); 
//            }
            
//            
          //Chen: Construct a diagonal matrix
            
            
            boolean Using1asDiagonal=true;
            double inpool_site_total_freq=0.0;
            double inpool_site_total_count=0.0;
            for (int l = 0; l < num_site_regional; l++) {
            	for (int i=0; i<inpool_site_freqs[l].length; i++ ) {
            		inpool_site_total_freq+= inpool_site_freqs[l][i];
            		inpool_site_total_count+=1.0;
            	}
            }
            if ((inpool_site_total_freq / inpool_site_total_count) <0.5) {
            	Using1asDiagonal=false;
            }
            System.out.println((inpool_site_total_freq) /inpool_site_total_count);
            
          //Chen: Construct a diagonal matrix
            
            for (int i=0;i< num_site_regional;i++) {
//            	if (Using1asDiagonal) {
//	            	String diagonal_hap =  diagonal_str(i, num_site_regional);
//	            	if (!haps_set.contains(diagonal_hap )) {
//	            		linked_haps.add(diagonal_hap);
//	            	}
//            	}
//            	if (!Using1asDiagonal) {
//	            	String diagonal_hap =  diagonal_0_str(i, num_site_regional);
//	            	if (!haps_set.contains(diagonal_hap )) {
//	            		linked_haps.add(diagonal_hap);
//	            	}
//            	}
            }
            
            int haps_2n = linked_haps.size();
            String[][] global_haps_string = new String[haps_2n][num_site_regional];
            String[] hap_IDs = new String[haps_2n];
            for (int h = 0; h < haps_2n; h++) {
                String curr_ID = "";
                String vc_str = linked_haps.get(h);
                String[] vc_arr  = vc_str.split("");
                // the length of vc_str may not reach num_site_regional; so put zeros in.
                for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                    global_haps_string[h][locus] = "0";
                    curr_ID += "0";
                }
                for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
                    global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
                    curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
                }
                hap_IDs[h] = curr_ID;
            }
//            for (int i=0;i<global_haps_string.length;i++ ) {
//            	for (int j=0;j <global_haps_string[i].length;j++ ) {
//            		System.out.print(global_haps_string[i][j]); 
//            	}
//            	System.out.println();
//            }
            //Chen: Construct a diagonal matrix
            
            double[] global_haps_freq = new double[haps_2n];
            Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
            return new HapConfig(
                global_haps_string,
                global_haps_freq,
                null,
                inpool_site_freqs,
                locusInfo,
                this.num_pools,
                hap_IDs,
                pool_IDs,
                this.dp.est_ind_pool,
                Using1asDiagonal);
        }
    
    
    
    public HapConfig generate_hapconfig_level_VIII_fail(
    		int[][] the_region, int region_index, String[] pool_IDs, int flag_index) 
    				throws IOException {

            int region_start = the_region[region_index][0];
            int region_end = the_region[region_index][1];
            int num_site_regional = region_end-region_start + 1;
            double[][] inpool_site_freqs = new double[num_site_regional][];
            LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];

            for (int l = 0; l < num_site_regional; l++) {
                locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
                inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
            }
            HashSet<String > haps_set =new HashSet<String >();
            ArrayList<String >  linked_haps= new ArrayList<String >();
            if ((this.num_regions_level_V%2==0) && (region_index==this.num_regions_level_VIII-1)){
            	for (int i =0;i< this.level_V_config[ this.level_V_config.length-1].num_global_hap;
            			i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_V_config[ this.level_V_config.length-1].
            				global_haps_string[i].length;j++) {
            			tmp=tmp+ this.level_V_config[ this.level_V_config.length-1].
                				global_haps_string[i][j];
            		}
            		linked_haps.add(tmp);
            	}
            }else {
            	ArrayList<String >  String_X= new ArrayList<String >();
            	ArrayList<String >  String_Y= new ArrayList<String >();
            	ArrayList<String >  String_Z= new ArrayList<String >();
            	int overlap= 1+ this.regions_level_V[2*region_index+1 ][1]- 
            			this.regions_level_VI[2*region_index +1][0];
//            	
            	int mismatch_cutoff=0;
            	for (int i =0;i< this.level_V_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_V_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_V_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_X.add(tmp);
            	}
            	for (int i =0;i< this.level_V_config[2*region_index+2].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_V_config[2*region_index+2].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_V_config[2*region_index+2].global_haps_string[i][j];
            		}
            		String_Y.add(tmp);
            	}
            	
            	for (int i =0;i< this.level_VI_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_VI_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_VI_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_Z.add(tmp);
            	}
//            	
            	for (int i =0;i< String_X.size();i++) {
            		for (int j =0;j< String_Y.size();j++) {
            			if (flag_index==-1) {
            				linked_haps.add(String_X.get(i)+String_Y.get(j) );
        					haps_set.add(String_X.get(i)+String_Y.get(j));
            			}else  if (flag_index!=-2){
            				if (!haps_set.contains(String_X.get(i)+String_Y.get(j)) ) {
            					haps_set.add(String_X.get(i)+String_Y.get(j));
            					linked_haps.add(String_X.get(i)+String_Y.get(j) );
            				}
            				for (int k =0;k< String_Z.size();k++) {
                				String[] comb_str = strmatch( String_X.get(i),
                						String_Y.get(j), String_Z.get(k),overlap,
                						flag_index);
                				for (int h=0;h<comb_str.length;h++ ) {
                					if (!haps_set.contains(comb_str[h])) {
                						linked_haps.add(comb_str[h]);
                						haps_set.add(comb_str[h]);
                					}
                				}
                			}
            			}else  if (flag_index==-2){
            				for (int k =0;k< String_Z.size();k++) {
	            				String[] comb_str = strmatch( String_X.get(i),
	            						String_Y.get(j), String_Z.get(k),overlap,
	            						this.level_VII_VIII_region_mismatch_tolerance);
	            				if (comb_str.length!=0) {
	            					String tmp_hap= String_X.get(i)+String_Y.get(j);
	            					if (!haps_set.contains(tmp_hap) ) {
	                					haps_set.add(tmp_hap);
	                					linked_haps.add(tmp_hap);
	                				}
	            					
		            				tmp_hap = tmp_hap.substring(0,this.regions_level_VI[2*region_index+1][0] -
		            						this.regions_level_V[2*region_index+1][0])
		            				+String_Z.get(k)
									+tmp_hap.substring(this.regions_level_VI[2*region_index+1][0]-
											this.regions_level_V[2*region_index+1][0]
											+ String_Z.get(k).length());
		            				if (!haps_set.contains(tmp_hap) ) {
	                					haps_set.add(tmp_hap);
	                					linked_haps.add(tmp_hap);
	                				}
	            				}
            				}
            			}
            		}
            	}
            }
            
            boolean Using1asDiagonal=true;
            double inpool_site_total_freq=0.0;
            double inpool_site_total_count=0.0;
            for (int l = 0; l < num_site_regional; l++) {
            	for (int i=0; i<inpool_site_freqs[l].length; i++ ) {
            		inpool_site_total_freq+= inpool_site_freqs[l][i];
            		inpool_site_total_count+=1.0;
            	}
            }
            if ((inpool_site_total_freq / inpool_site_total_count) <0.5) {
            	Using1asDiagonal=false;
            }
            System.out.println((inpool_site_total_freq) /inpool_site_total_count);
            
            
            for (int i=0;i< num_site_regional;i++) {
            	if (Using1asDiagonal) {
	            	String diagonal_hap =  diagonal_str(i, num_site_regional);
	            	if (!haps_set.contains(diagonal_hap )) {
	            		linked_haps.add(diagonal_hap);
	            	}
            	}
            	if (!Using1asDiagonal) {
	            	String diagonal_hap =  diagonal_0_str(i, num_site_regional);
	            	if (!haps_set.contains(diagonal_hap )) {
	            		linked_haps.add(diagonal_hap);
	            	}
            	}
            }
            
            int haps_2n = linked_haps.size();
            String[][] global_haps_string = new String[haps_2n][num_site_regional];
            String[] hap_IDs = new String[haps_2n];
            for (int h = 0; h < haps_2n; h++) {
                String curr_ID = "";
                String vc_str = linked_haps.get(h);
                String[] vc_arr  = vc_str.split("");
                for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                    global_haps_string[h][locus] = "0";
                    curr_ID += "0";
                }
                for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
                    global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
                    curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
                }
                hap_IDs[h] = curr_ID;
            }
            
            double[] global_haps_freq = new double[haps_2n];
            Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
            return new HapConfig(
                global_haps_string,
                global_haps_freq,
                null,
                inpool_site_freqs,
                locusInfo,
                this.num_pools,
                hap_IDs,
                pool_IDs,
                this.dp.est_ind_pool,
                Using1asDiagonal);
        }
    
    
    
    public HapConfig generate_hapconfig_level_VIII(
            int[][] the_region, int region_index, String[] pool_IDs) throws IOException {

            int region_start = the_region[region_index][0];
            int region_end = the_region[region_index][1];
            int num_site_regional = region_end-region_start + 1;
            double[][] inpool_site_freqs = new double[num_site_regional][];
            LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];

            for (int l = 0; l < num_site_regional; l++) {
                locusInfo[l] = this.in_pool_sites_freq_anno.loci_annotations[l + region_start];
                inpool_site_freqs[l] = this.in_pool_sites_freq_anno.inpool_freqs[l + region_start];
            }
            HashSet<String > haps_set =new HashSet<String >();
            ArrayList<String >  linked_haps= new ArrayList<String >();
            if ((this.num_regions_level_V%2==0) && (region_index==this.num_regions_level_VIII-1)){
            	for (int i =0;i< this.level_V_config[ this.level_V_config.length-1].num_global_hap;
            			i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_V_config[ this.level_V_config.length-1].
            				global_haps_string[i].length;j++) {
            			tmp=tmp+ this.level_V_config[ this.level_V_config.length-1].
                				global_haps_string[i][j];
            		}
            		linked_haps.add(tmp);
            	}
            }else {
            	ArrayList<String >  String_X= new ArrayList<String >();
            	ArrayList<String >  String_Y= new ArrayList<String >();
            	ArrayList<String >  String_Z= new ArrayList<String >();
            	int overlap= 1+ this.regions_level_V[2*region_index+1 ][1]- 
            			this.regions_level_VI[2*region_index +1][0];
//            	
            	int mismatch_cutoff=0;
            	for (int i =0;i< this.level_V_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_V_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_V_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_X.add(tmp);
            	}
            	for (int i =0;i< this.level_V_config[2*region_index+2].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_V_config[2*region_index+2].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_V_config[2*region_index+2].global_haps_string[i][j];
            		}
            		String_Y.add(tmp);
            	}
            	
            	for (int i =0;i< this.level_VI_config[2*region_index+1].num_global_hap;i++) {
            		String tmp ="";
            		for (int j =0; j< this.level_VI_config[2*region_index+1].global_haps_string[i].length;
            				j++) {
            			tmp=tmp+ this.level_VI_config[2*region_index+1].global_haps_string[i][j];
            		}
            		String_Z.add(tmp);
            	}
//            	
            	for (int i =0;i< String_X.size();i++) {
            		for (int j =0;j< String_Y.size();j++) {
            			for (int k =0;k< String_Z.size();k++) {
            				String[] comb_str = strmatch( String_X.get(i),
            						String_Y.get(j), String_Z.get(k),overlap,
            						this.level_VII_VIII_region_mismatch_tolerance);
            				if   ((comb_str.length>0) &&
            						(!haps_set.contains(String_X.get(i)+  String_Y.get(j)))) {
            					linked_haps.add(String_X.get(i)+  String_Y.get(j));
        						haps_set.add(String_X.get(i)+  String_Y.get(j));
        						break;
            				}
            				for (int h=0;h<comb_str.length;h++ ) {
            					if (!haps_set.contains(comb_str[h])) {
            						linked_haps.add(comb_str[h]);
            						haps_set.add(comb_str[h]);
            					}
            				}
            			}
            		}
            	}
            }
            
            boolean Using1asDiagonal=true;
            double inpool_site_total_freq=0.0;
            double inpool_site_total_count=0.0;
            for (int l = 0; l < num_site_regional; l++) {
            	for (int i=0; i<inpool_site_freqs[l].length; i++ ) {
            		inpool_site_total_freq+= inpool_site_freqs[l][i];
            		inpool_site_total_count+=1.0;
            	}
            }
            if ((inpool_site_total_freq / inpool_site_total_count) <0.5) {
            	Using1asDiagonal=false;
            }
            System.out.println((inpool_site_total_freq) /inpool_site_total_count);
            
          //Chen: Construct a diagonal matrix
            
            for (int i=0;i< num_site_regional;i++) {
//            	if (Using1asDiagonal) {
//	            	String diagonal_hap =  diagonal_str(i, num_site_regional);
//	            	if (!haps_set.contains(diagonal_hap )) {
//	            		linked_haps.add(diagonal_hap);
//	            	}
//            	}
//            	if (!Using1asDiagonal) {
//	            	String diagonal_hap =  diagonal_0_str(i, num_site_regional);
//	            	if (!haps_set.contains(diagonal_hap )) {
//	            		linked_haps.add(diagonal_hap);
//	            	}
//            	}
            }
            
            int haps_2n = linked_haps.size();
            String[][] global_haps_string = new String[haps_2n][num_site_regional];
            String[] hap_IDs = new String[haps_2n];
            for (int h = 0; h < haps_2n; h++) {
                String curr_ID = "";
                String vc_str = linked_haps.get(h);
                String[] vc_arr  = vc_str.split("");
                // the length of vc_str may not reach num_site_regional; so put zeros in.
                for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                    global_haps_string[h][locus] = "0";
                    curr_ID += "0";
                }
                for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
                    global_haps_string[h][locus] = vc_arr[locus - num_site_regional + vc_arr.length];
                    curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
                }
                hap_IDs[h] = curr_ID;
            }
//            for (int i=0;i<global_haps_string.length;i++ ) {
//            	for (int j=0;j <global_haps_string[i].length;j++ ) {
//            		System.out.print(global_haps_string[i][j]); 
//            	}
//            	System.out.println();
//            }
            //Chen: Construct a diagonal matrix
            
            double[] global_haps_freq = new double[haps_2n];
            Arrays.fill(global_haps_freq, 1.0 / (double) haps_2n);
            return new HapConfig(
                global_haps_string,
                global_haps_freq,
                null,
                inpool_site_freqs,
                locusInfo,
                this.num_pools,
                hap_IDs,
                pool_IDs,
                this.dp.est_ind_pool,
                Using1asDiagonal);
        }
    
    
    /**
     *Chen: Using the global_haps from possible paths instead of gc 
     *  @param region
     *  @param level
     *  @param r_index
     *  @param dir_prefix intermediate folder
     *  @return
     *  @throws FileNotFoundException
     */
    public HapConfig generate_hapconfig_gc_regional_lasso( 
        String[] pool_IDs,
        String[] vef_files,
        int[] region,
        int level,
        int r_index,
        String aem_fail_lasso_path)  // i 
            throws FileNotFoundException {

        new File(aem_fail_lasso_path).mkdir();
        int start = region[0]; 
        int end = region[1];
        int num_site_regional = end - start + 1;
        // collect all regional haps from raw GC output. 
        
        HashMap<String, Double> hap_tracker = new HashMap<String, Double>();
        
        int haps_2n = (int) Math.pow(2, num_site_regional);
        String[][] global_haps_string_tmp  = new String[haps_2n][num_site_regional];
        String[] hap_IDs_tmp = new String[haps_2n];
        double[] global_haps_freq_tmp  = new double [haps_2n];
        double total_freq =0;
        for (int h = 0; h < haps_2n; h++) {
        	String curr_ID = "";
        	String vc_str = Integer.toBinaryString(h);
        	String[] vc_arr  = vc_str.split("");
        	// the length of vc_str may not reach num_site_regional; so put zeros in.
        	for (int locus = 0; locus < num_site_regional - vc_arr.length; locus++) {
                global_haps_string_tmp[h][locus] = "0";
                curr_ID += "0";
            }
        	for (int locus = num_site_regional - vc_arr.length; locus<num_site_regional; locus++) {
        		global_haps_string_tmp[h][locus] = vc_arr[locus - num_site_regional 
        		                                          + vc_arr.length];
        		curr_ID += vc_arr[locus - num_site_regional + vc_arr.length];
        	}
        	hap_IDs_tmp[h] = curr_ID;
        	for (int p=0; p< this.num_pools;p++) {
        		double freq = freq_cal(this.loci_link_freq[p],curr_ID, start,end );
        		global_haps_freq_tmp[h]+= freq ;
        		if (freq >  1 / (double) haps_2n){
        			total_freq+= freq;
        		}
        	}
        }
        for (int h = 0; h < global_haps_freq_tmp.length; h++) {
        	if (global_haps_freq_tmp[h] > 1 / (double) haps_2n) {
        		String tmp= "";
        		for (int j=0; j < global_haps_string_tmp[h].length;j++) {
            		tmp =tmp+ global_haps_string_tmp[h][j];
           		}
        		hap_tracker.put(tmp, global_haps_freq_tmp[h]/ total_freq ); 
        	}
        }
        
//        for (int h = 0; h < this.num_haps_gc; h++) {
//            String curr_vc = String.join("",Arrays.copyOfRange(this.global_haps_gc[h],start,end+1));
//            if (!hap_tracker.containsKey(curr_vc)) {
//                hap_tracker.put(curr_vc, this.global_gc_freq[h]); 
//            } else {
//                double new_freq = hap_tracker.get(curr_vc) + this.global_gc_freq[h];
//                hap_tracker.put(curr_vc, new_freq);
//            }
//        }
        
        int num_hap_regional = hap_tracker.size();
        String[][] reg_haps_string = new String[num_hap_regional][num_site_regional];
        int hap_index = 0;
        String[] hap_IDs = new String[num_hap_regional];
        double[] global_haps_freq = new double[num_hap_regional];
        for (String curr_vc : hap_tracker.keySet()) {
            String[] tmp = curr_vc.split("");
            for (int s = 0; s < num_site_regional; s++) reg_haps_string[hap_index][s] = tmp[s];
            global_haps_freq[hap_index] = hap_tracker.get(curr_vc) / num_pools;
            hap_IDs[hap_index] = "h"+Integer.toString(hap_index);
            hap_index++;
        }

        double[][] inpool_site_freqs = new double[num_site_regional][];
        LocusAnnotation[] locusInfo = new LocusAnnotation[num_site_regional];
        for (int s = 0; s < num_site_regional; s++) {
            locusInfo[s] = this.in_pool_sites_freq_anno.loci_annotations[s + start];
            inpool_site_freqs[s] = this.in_pool_sites_freq_anno.inpool_freqs[s + start];
        }

        HapConfig gc_regional_haps = new HapConfig(
            reg_haps_string,
            global_haps_freq,
            null,
            inpool_site_freqs,
            locusInfo,
            this.num_pools,
            hap_IDs,
            pool_IDs,
            this.dp.est_ind_pool);

        HapLASSO regional_lasso = new HapLASSO( 
            -1,    // Specified by pool_index==-1, this is the *multi_pool* implementation of LASSO.
            null,  // pool_IDs is null as it is the multi_pool LASSO.
            dp.lasso_regional_lambda,
            gc_regional_haps,
            0,
            dp.lasso_regional_memory,
            aem_fail_lasso_path + "/"+dp.project_name+"_level_" + level + "_region_" + r_index + ".lasso_in", 
            1, 1 );

        regional_lasso.estimate_frequencies_lasso(null, vef_files, dp.lasso_weights);

        HapConfig final_reg_haps = regional_lasso.hapOut(pool_IDs);

        // filter out low-frequent haps 
        boolean[] list_remove_haps = new boolean[final_reg_haps.num_global_hap];      
        int num_remove_hap = 0;
        for (int h = 0; h < final_reg_haps.num_global_hap; h++) {
            if (final_reg_haps.global_haps_freq[h] < dp.lasso_regional_cross_pool_cutoff) {
                list_remove_haps[h] = true;
                num_remove_hap++;
            }
        }
        double actual_cutoff = dp.lasso_regional_cross_pool_cutoff;
        
        // If too many of them are below the regional frequency minimum,
        // i.e., too few haps are remained 
        if (dp.lasso_hapset_size_min > final_reg_haps.num_global_hap - num_remove_hap) {
            if(final_reg_haps.num_global_hap >= dp.lasso_hapset_size_min) {
                list_remove_haps = Algebra.permute_sort_and_remove(
                    final_reg_haps.global_haps_freq.clone(), 
                    dp.lasso_hapset_size_min);
            }            
        }
        // Or, if too many haps are remained...
        else if (dp.lasso_hapset_size_max < final_reg_haps.num_global_hap - num_remove_hap) {
            list_remove_haps = Algebra.permute_sort_and_remove(
                final_reg_haps.global_haps_freq.clone(), 
                dp.lasso_hapset_size_max);
        }
        num_remove_hap = 0;
        actual_cutoff=1;
        for(int h=0; h < final_reg_haps.num_global_hap; h++) {
            if(list_remove_haps[h]) {
                num_remove_hap++;
            } else {
                //this.initial_Haps.global_haps_freq[h] = freq[h];
                if(actual_cutoff > final_reg_haps.global_haps_freq[h]) { 
                    actual_cutoff = final_reg_haps.global_haps_freq[h];
                }
            }
        }     
        final_reg_haps.remHaps(list_remove_haps, num_remove_hap);
        System.out.println("Of the "
            + final_reg_haps.num_global_hap
            + " regional haplotypes, "
            + num_remove_hap
            + " were removed. "
            + "The frequency cutoff was " + actual_cutoff + ".");        
//        if (Double.isNaN(actual_cutoff)) {
//            double[] tmp = new double[final_reg_haps.num_global_hap];
//            for (int h = 0; h < final_reg_haps.num_global_hap; h++) {
//                tmp[h] = 1.0 / ((double) final_reg_haps.num_global_hap);
//            }
//            final_reg_haps.global_haps_freq = tmp;
//        }
        return final_reg_haps;
    }
}
