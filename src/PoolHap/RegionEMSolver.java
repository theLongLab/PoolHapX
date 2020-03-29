package PoolHap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import PoolHap.Parameters;

public class RegionEMSolver {
    // Redundant variables for convenience.
    public int num_loci; // number of SNPs in this region
    public int num_pools; // number of pools under study

    // Main structures.
    public HapConfig initial_Haps;
    public boolean failure = false;
    public HapConfig final_Haps;
    public HashSet<Integer> error_index= new HashSet<Integer>() ;
    
    public double [] aem_freq;
    public String [][] aem_haps;
    public double aem_min_diff;
    public int aem_level;
    public double fixed_diff;
    public double hc_cufoff;
    public Parameters gp;
    public int region;
    
    public double if_0_0;
    public double if_Denominator_0;
    
    
    

    /**
     *  The constructor for the PoolSolver object i.e.: the core algorithm.
     *
     *  @param Estimated variant frequencies, positions of the primitive loci, list of each pool's
     *  VEFs.
     *
     *  @return To set up and adjust variant composition and inter/intra-pool frequencies of
     *  curr_Haps.
     */
    public RegionEMSolver(HapConfig hap_config, String parameter_file , int level, int reg, 
    		double diff) throws Exception {
//    	this.initial_Haps.
        this.initial_Haps = hap_config;
        this.num_loci = this.initial_Haps.num_loci;
        this.num_pools = this.initial_Haps.num_pools;
        
        this.initial_Haps.regional_EM_solver = this;
        this.aem_level = level;
        this.fixed_diff = diff;
        this.region = reg;
        analyze_a_region_aem(parameter_file);
        
      
        
        // analyze_a_region_mcmc(parameter_file);
    }

    /**
     *  The main method for running EM on haplotypes to estimate their global frequency. Inspired by
     *  the AEM function from Kuk et al., 2009.
     *
     *  @param this.initial_Haps, the set of possible sub-haplotypes in a specific region.
     *  aem_parameters include est_ind_pool, epsilon, rare_cutoff.
     *
     *  @return Updated guess of the overall (between-pool frequencies of the sub-haplotypes)
     *  frequencies i.e.: updated this.initial_Haps object.
     */
    
    public int num_of_alternate (String [] hap ) throws IOException {
    	int count= 0;
    	for (int i =0; i< hap.length;i++) {
    		if (hap[i].equals("1")) {
    			 count++;
    		}
    	}
    	return count;
    }
    
    public void track_gold_standard (String [][] hap, double [] hap_freq ) 
    		throws IOException {
    	String []compare_haps =  new String [hap_freq.length ];
    	double [] compare_freq =  new double [hap_freq.length ];
    	for (int i =0;i< hap_freq.length;i++) {
    		compare_freq[i] = hap_freq[i];
    		String tmp ="";
    		for (int j =0; j<hap[i].length;j++ ) {
    			tmp=tmp+ hap[i][j];
    		}
    		compare_haps[i] =tmp;
    	}
    	for (int i =0;i< compare_freq.length;i++) {
    		for (int j=i;j< compare_freq.length;j++) {
    			if (compare_freq[i]< compare_freq[j]) {
    				double tmp_freq= compare_freq[i];
    				compare_freq[i] =compare_freq[j];
    				compare_freq[j]= tmp_freq;
    				String tmp_hap = compare_haps[i];
    				compare_haps[i]=compare_haps[j];
    				compare_haps[j]= tmp_hap;
    			}
    		}
    	}
    	
    	
    	
    	String dc_file = gp.inter_dir + gp.project_name + "_dc_plan.txt";
    	
    	ArrayList<Integer > level_I_start = new ArrayList<Integer>();
		ArrayList<Integer > level_I_end = new ArrayList<Integer>();
		ArrayList<Integer > level_II_start = new ArrayList<Integer>();
		ArrayList<Integer > level_II_end = new ArrayList<Integer>();
		int count = 0;
		BufferedReader bufferedreader = new BufferedReader(new FileReader(dc_file));
        String line = "";
        while ((line = bufferedreader.readLine()) != null) {
        	line =line.replace("\n", "").replace("\r", "");
//        	System.out.println(line);
        	if (count== 5) {
        		String[] tmp = line.split("\t");
        		for (int i = 0; i < tmp.length; i++) {
        			String[] tmp2= tmp[i].split(":");
        			level_I_start.add(Integer.parseInt(tmp2[0] )) ;
        			level_I_end.add(Integer.parseInt(tmp2[1] )) ;
        		}
        	}else if (count==7) {
        		String[] tmp = line.split("\t");
        		for (int i = 0; i < tmp.length; i++) {
        			String[] tmp2= tmp[i].split(":");
        			level_II_start.add(Integer.parseInt(tmp2[0] )) ;
        			level_II_end.add(Integer.parseInt(tmp2[1] )) ;
        		}
        	}
        	count ++;
        }
        bufferedreader.close();
        
        ArrayList<String > gold_haps =new ArrayList<String > ();
        String gold_file = this.gp.gold_dir+"/"+this.gp.project_name+"_haps.inter_freq_vars.txt";
        BufferedReader bufferedreader2 = new BufferedReader(new FileReader(gold_file));
        line = "";
        ArrayList<ArrayList<String > > tmp_2d = new ArrayList<ArrayList<String>>();
        while ((line = bufferedreader2.readLine()) != null) {
        	line =line.replace("\n", "").replace("\r", "");
        	if ((!line.startsWith("Hap_ID")) &&  (!line.startsWith("Freq") )){
        		String[] tmp = line.split("\t");
        		ArrayList<String> tmp_arr = new ArrayList<String>();
        		for (int j = 1; j < tmp.length; j++) {
        			tmp_arr.add(tmp[j]);
        		}
        		tmp_2d.add(tmp_arr);
        	}
        }
        bufferedreader2.close();
        for (int j = 0; j < tmp_2d.get(0).size(); j++) {
        	String tmp_str="";
        	for (int k = 0; k < tmp_2d.size(); k++) {
        		tmp_str=tmp_str+ tmp_2d.get(k).get(j);
        	}
        	gold_haps.add(tmp_str);
        }
        
        

        int total_mismatches=0;
    	for (int j = 0; j < gold_haps.size(); j++) {
        	int  max_mismatch=  100;
        	int index =0;
        	for (int k = 0; k < compare_haps.length; k++) {
        		if( NumofMismatch (gold_haps.get(j).substring(level_I_start.get(this.region), 
        				level_I_end.get(this.region)+1), compare_haps[k] ) <  max_mismatch) {
        			max_mismatch= NumofMismatch (gold_haps.get(j).substring(level_I_start.get(this.region), 
            				level_I_end.get(this.region)+1), compare_haps[k] );
        			index =k;
        		}
        	}
        	System.out.println(gold_haps.get(j).substring(level_I_start.get(this.region), level_I_end.get(this.region)+1)
        			+"\tFor Level III: Region "+ Integer.toString(this.region)+": Mismatch: "
        			+Integer.toString(max_mismatch)+"\t"+ compare_haps[index] +"\t"+index+"\t"+
        			compare_freq[index]);
        	total_mismatches  += max_mismatch;
        }
    	
    	
    	total_mismatches=0;
    	for (int j = 0; j < gold_haps.size(); j++) {
        	int  max_mismatch=  100;
        	int index =0;
        	for (int k = 0; k < 50; k++) {
        		if( NumofMismatch (gold_haps.get(j).substring(level_I_start.get(this.region), 
        				level_I_end.get(this.region)+1), compare_haps[k] ) <  max_mismatch) {
        			max_mismatch= NumofMismatch (gold_haps.get(j).substring(level_I_start.get(this.region), 
            				level_I_end.get(this.region)+1), compare_haps[k] );
        			index =k;
        		}
        	}
//        	System.out.println(gold_haps.get(j).substring(level_I_start.get(this.region), level_I_end.get(this.region)+1)
//        			+"\tFor Level III: Region "+ Integer.toString(this.region)+": Mismatch: "
//        			+Integer.toString(max_mismatch)+"\t"+ compare_haps[index] +"\t"+index+"\t"+
//        			compare_freq[index]);
        	total_mismatches  += max_mismatch;
        }
    	
    	
    	
    	System.out.println("************************\t"+ total_mismatches);

    }
    
    
    public int NumofMismatch (String x, String y) throws IOException {
		int mismatch= 0;
		for (int i=0;i < x.length();i++){
			if (!x.substring(i, i+1).equals( y.substring(i, i+1))) {
				mismatch++ ;
			}
		}
		return mismatch;
	}
    
    public boolean check_hap (String [][] hap, double [] hap_freq , boolean flag) 
    		throws IOException {
    	
    	double diagonal_matrix_freq =0.0;
    	double total_freq=0.0001;
    	if (hap[0].length<5) {
    		return true;
    	}
    	for (int i=0;i< hap.length;i++) {
    		if (( flag) && (num_of_alternate(hap[i])==1 )) {
    			diagonal_matrix_freq+= hap_freq[i];
    		} else if (( !flag) && (num_of_alternate(hap[i])==(hap[i].length-1 ) )) {
    			diagonal_matrix_freq+= hap_freq[i];
    		}
    		total_freq+= hap_freq[i];
    		
    	}
    	if (( diagonal_matrix_freq/ total_freq)> 0.9) {
    		return false;
    	}
    	
    	HashSet<Double> freq_dict = new HashSet<Double>  ();
    	
    	double min_freq = 0.00000001;
    	int exceed_min_freq_count=0;
    	
    	
    	for (int i=0;i< hap.length;i++) {

    		if (( !flag) || (num_of_alternate(hap[i])!=1 )) {
    			freq_dict.add( hap_freq[i]);
    			if  (hap_freq[i]> min_freq) {
    				exceed_min_freq_count++;
    			}
    		}
    		if (( flag) || (num_of_alternate(hap[i])!=(hap[i].length-1 )  )) {
    			freq_dict.add( hap_freq[i]);
    			if  (hap_freq[i]> min_freq) {
    				exceed_min_freq_count++;
    			}
    		}
    		
    	}
    	if ((exceed_min_freq_count<5) || (freq_dict.size()<5)) {
    		return false;
    	}
    	
    	return true;
    }
    
    public double freq_dist (double [] x, double [] y) {
    	double n =(double) x.length;
    	double value=0;
    	for (int i=0;i< x.length;i++) {
    		value += Math.abs(x[i]-y[i]);
    	}
    	return value/ n;
    }
    
    public double [] var_freq_cal (String [][] hap, double [] hap_freq ) throws IOException {
    	double [] var_freq =new double [hap[0].length ];
    	for (int i=0;i< var_freq.length; i++) {
    		double value =0;
    		for (int j=0; j< hap.length; j++) {
    			if (hap[j][i].equals("1")) {
    				value =value+ hap_freq[j];
    			}
    		}
    		var_freq[i]=value;
    	}
    	return var_freq;
    }
    
    public double [] var_freq_cal2 ( double [][] site_freq ) throws IOException {
    	double [] var_freq =new double [site_freq.length ];
    	for (int i=0;i< var_freq.length; i++) {
    		double value =0;
    		for (int j=0; j< site_freq[i].length; j++) {
    			value =value+ site_freq[i][j];
    		}
    		var_freq[i]=value/ (double)site_freq[i].length ;
    	}
    	return var_freq;
    }
    
    public double [] outlier_rh( double [] rh ) throws IOException {
    	double [] rh_new  =rh.clone();
    	
    	for (int i=0;i< rh.length; i++) {
    		for (int j=i; j< rh.length; j++) {
    			if (rh[j]> rh[i]) {
    				double value =0;
    				value= rh[i];
    				rh[i]=rh[j];
    				rh[j]= value;
    			}
    		}
    	}
    	double max_value=0.0;
    	boolean flag= false;
    	for (int i=0;i< (rh.length-1); i++) {
    		if ((rh[i]> (10.0* rh[i+1]))  && (i<this.num_loci)){
    			max_value = rh[i+1];
    			flag=true;
    		}
    	}
    	
    	for (int i=0;i< rh.length; i++) {
    		if ((rh_new[i]> max_value)  && (flag)){
    			rh_new[i] = max_value;
    		}
    	}
    	return rh_new;
    }
    
    
    public void analyze_a_region_aem(String parameter_file) throws IOException {	
    	this.aem_min_diff =1.0;
        Parameters aem_parameters = new Parameters(parameter_file);
        this.gp = aem_parameters;
        
        this.hc_cufoff= aem_parameters.hc_similarity_cutoff;
        // The current estimate of global haplotype frequencies.
        double[] freq = this.initial_Haps.global_haps_freq.clone();
        
//        System.out.println("Number of potential haplotypes\t"+ freq.length);
//        for (int i = 0; i  < freq.length; i++) {    		
//        	System.out.print(freq[i]  );
//        	System.out.print("\t" );
//        	for (int j = 0; j  < this.initial_Haps.global_haps_string[i].length; j++) { 
//        		System.out.print(this.initial_Haps.global_haps_string[i][j]  );
//        	}
//        	System.out.print("\n" );
//        }
        // The importance factor of each haplotype is the previous iteration's estimate of its
        // global frequency.
        double[] Rh = freq.clone();
        double[] freq_bk = freq.clone();
        this.initial_Haps.update_sigma_mu_logL();
        int  iter_coefficient = 10;
//        if (this.num_loci<20) {
//        	 iter_coefficient =6;
//        }else  if (this.num_loci> 20) {
//        	iter_coefficient =100;
//        }
        int max_iter =iter_coefficient * this.num_loci;
        double min_freq_cutoff= 1.0;
        int iter=0;
        int diagonal_matrix_freq_index =0;
        double iter_aem_min_diff=1.0;
//        double[] diagonal_matrix_freq =new double [] {0.0002,0.0, -0.0001};
        double[] diagonal_matrix_freq =new double [] {0.0002};
      
//        max_iter= 50;
        while ((iter< max_iter)  &&
        		(diagonal_matrix_freq_index < diagonal_matrix_freq.length ) )  {
//        	System.out.println("______________________");
//        	System.out.println(iter+"\t"+"______________________");
        	double total_sigma =0;
        	for (int j = 0; j < this.initial_Haps.sigma.length; j++) {
        		for (int k = 0; k < this.initial_Haps.sigma[j].length; k++) {
        			total_sigma+= this.initial_Haps.sigma[j][k];
        		}
        	}
        	
        	
        	
        	
        	double total_sigma_cufoff= 0.02* (double)this.num_loci* (double)this.num_loci;
        	if (total_sigma< total_sigma_cufoff) {
        		double multiple = total_sigma_cufoff/ total_sigma;
        		for (int j = 0; j < this.initial_Haps.sigma.length; j++) {
            		for (int k = 0; k < this.initial_Haps.sigma[j].length; k++) {
            			this.initial_Haps.sigma[j][k] = this.initial_Haps.sigma[j][k]* multiple;
            		}
            	}
        	}
        	
        	
            SingularValueDecomposition svd = new SingularValueDecomposition(
                MatrixUtils.createRealMatrix(this.initial_Haps.sigma));

            RealMatrix svd_inv = svd.getSolver().getInverse();
            
            // Step 1b) Calculate i) by taking the distance of alternate alleles from the average
            // haplotype of each pool.
            // TODO: [ReconEP]:: extract to helper?
            double[][] dist_to_mu = new double[this.num_pools][this.num_loci];
            // dist_to_mu = sweep(A / (2*N), 2, omega, "-");
            for (int u = 0; u < this.num_pools; u++) {
                for (int v = 0; v < this.num_loci; v++) {
                    dist_to_mu[u][v] = this.initial_Haps.inpool_site_freqs[v][u]
                        - this.initial_Haps.mu[v];
                }
            } 
            
//            for (int u = 0; u < this.num_pools; u++) {
//                for (int v = 0; v < this.num_loci; v++) {
//                   System.out.print(dist_to_mu[u][v] );;
//                }
//                System.out.println();
//            } 
//            
            
//            BufferedWriter bw = new BufferedWriter(new FileWriter("111.txt", true));

            
            double[] IF = new double[this.initial_Haps.num_global_hap];
            for (int j = 0; j < this.initial_Haps.num_global_hap; j++) {
                double[] hh = this.initial_Haps.global_haps[j];
                // rh1 = exp(-1 / (4 * N) * t(omega - h) %*% svd_inv %*% (omega - h))
                // The approximation of the importance factor of each haplotype, according to
                // Appendix B of Zhang et al. 2008.
                
                double rh1 = Math.exp(-1 / (2.0 * aem_parameters.est_ind_pool)
                    * Algebra.quadratic_form(
                        Algebra.minus(this.initial_Haps.mu, hh),
                        svd_inv.getData()));
//                System.out.println(rh1);
                // QUESTION: is this.aem_parameters.est_ind_pool making rh1 a lot smaller than it
                // has to be?
                RealMatrix dist_mtx = MatrixUtils.createRealMatrix(dist_to_mu);
                // rh2 = exp(-dist_to_mu %*% svd_inv %*% (omega-h)
                //     - diag(dist_to_mu %*% svd_inv %*% t(dist_to_mu) / 2))

                
                double[] rh2 = Algebra.exp(
                    Algebra.minus(svd_inv.preMultiply(dist_mtx)
                        .operate(Algebra.minus(hh, this.initial_Haps.mu)),
                    Algebra.diag(dist_mtx.multiply(svd_inv)
                        .multiply(dist_mtx.transpose())
                        .scalarMultiply(0.5)
                        .getData())));
                
                
//                rh2= outlier_rh(rh2);
//                Algebra.normalize_ditribution(rh2);
                double rh_min = this.gp.if_0_0;
                double rh_max = this.gp.if_Denominator_0 ;
                double rh2_mean = Algebra.mean(rh2);
                if ( (Double.isInfinite( rh2_mean))   
                		||  (rh2_mean> rh_max)  ){
                	rh2_mean = rh_max;
                }else if (( Double.isNaN( rh2_mean))   || (rh2_mean< rh_min)) {
                	rh2_mean = rh_min;
                }
                
                if ( (Double.isInfinite( rh1) ) ||  (rh1>rh_max )     ){
                	rh1 = rh_max;
                }else if (( Double.isNaN( rh1))   || (rh1< rh_min)) {
                	rh1 = rh_min;
                }
                
//                if (j==0){
//                	System.out.print("\nFreq:");
//                	if (freq.length>50) {
//                	for (int i=0;i< 50; i++) {
//                		System.out.print(freq[i]+"\t");
//    	            }
//                	}
//	                System.out.print("\nrh2:");
//	                for (int i=0;i< rh2.length; i++) {
//	                	System.out.println(this.num_loci+"\t"+rh1+"\t"+rh2[i]);
//	                }
//                }	
                
                
                double rh = rh1 * rh2_mean;
//                if (this.num_loci>20) {
//                	bw.write("\nrh\t"+rh1+"\t"+Algebra.mean(rh2));
//                }
//                rh = rh1 ;
                IF[j] = rh;
            }
            
            double delta = Algebra.sum(Algebra.abs(Algebra.minus(Rh, IF))); // sum(abs(Rh - IF)
            
            // sum(abs(IF / Rh - 1))
            double delta1 = Algebra.sum(Algebra.abs(Algebra.add(Algebra.divide(Rh,  IF), -1)));
            
            double sum = Algebra.sum(IF);// debug by Chen 2019-10
            double freq_inpool=0.0;
            double num_freq=0.000001;
            
//            for (int i=0;i< IF.length; i++) {
//            	if (IF[i]>100.0) {
//            		IF[i]= 100.0;
//            	}
//            	if ((Double.isNaN(IF[i]))  || (IF[i]<1e-10))  {
//            		IF[i]= 1e-10;
//            	}
//            }
            
//            for (int i=0;i< IF.length; i++) {
//            	if (IF[i]>0) {
//            		num_freq=num_freq+1;
//            	}
//            	if (IF[i]>    (sum*0.1)) {
//            		IF[i]=sum*0.1 ;
//            		freq_inpool= freq_inpool+IF[i]*0.9;
//            	}
//            }
//            for (int i=0;i< IF.length; i++) {
//            	if (IF[i]>0) {
//            		IF[i]= IF[i]+ freq_inpool/ num_freq;
//            	}
//            }
//            
            
            Rh = IF;
//            BufferedWriter bw = new BufferedWriter(new FileWriter("/export/qlong/chencao/Work/poolhapx/sim/111.txt", true));
//            if (freq.length>40) {
//            	bw.write("\niter:"+ iter+"\t"+IF.length +"\n");
//            	bw.write("\nfreq_ini:\n");
//	            for (int i=0;i< freq.length; i++) {
//	            	bw.write(freq[i]+"\t");
//	            }
//                       
//	            bw.write("\nIF:\n");
//	            for (int i=0;i< IF.length; i++) {
//	            	bw.write(IF[i]+"\t");
//	            }
//            }
            
            double[] p_new = Algebra.times(Rh, freq);
//            if (freq.length>40) {
//            	bw.write("\np_new:\n");
//	            for (int i=0;i< p_new.length; i++) {
//	            	bw.write(p_new[i]+"\t");
//	            }
//            }

            Algebra.normalize_ditribution(p_new);
//            if (freq.length>40) {
//            	bw.write("\np_new_normalize:\n");
//	            for (int i=0;i< p_new.length; i++) {
//	            	bw.write(p_new[i]+"\t");
//	            }
//            }
            Algebra.rmlow_and_normalize(p_new, aem_parameters.aem_zero_cutoff);
            
            freq = p_new.clone();
            
//            if (freq.length>40) {
//            	bw.write("\nfreq_final:\n");
//	            for (int i=0;i< freq.length; i++) {
//	            	bw.write(freq[i]+"\t");
//	            }
//            }
//            bw.close();
            
//          if (freq.length>0) {
//        	System.out.print("\nfreq_final:\n");
//            for (int i=0;i< freq.length; i++) {
//            	System.out.print(freq[i]+"\t");
//            }
//          }
            
            
            
            for (double f : freq) {
                if (Double.isNaN(f))  {
                	freq= freq_bk.clone();
                	this.failure=false;
    	            this.initial_Haps.global_haps_freq = freq.clone();
    	            break;
                	
                }
            }
            
            double [] predicted_var_freq =new double [ this.initial_Haps.num_loci];
        	double [] real_var_freq =new double [ this.initial_Haps.num_loci];
        	predicted_var_freq= var_freq_cal(this.initial_Haps.global_haps_string,freq);
        	real_var_freq =var_freq_cal2(this.initial_Haps.inpool_site_freqs);
        	double freq_dist_value= freq_dist( predicted_var_freq,real_var_freq );
            if (freq_dist_value<min_freq_cutoff  ) {
            	System.out.println(min_freq_cutoff);
            	min_freq_cutoff = freq_dist_value;
            	freq_bk= freq.clone();
            }
            
           
            if (delta < aem_parameters.aem_epsilon || delta1 < aem_parameters.aem_epsilon) {
	            this.failure=false;
	            this.initial_Haps.global_haps_freq = freq.clone();
//	            break;
            }
            
            
            iter++;
         

            this.initial_Haps.global_haps_freq = freq.clone();
            
// Chen: Ensure a non-singular matrix   
            
            
            
            this.initial_Haps.update_sigma_mu_logL();
        }
//---------------------------------------- AEM end       
        
        
        
        
        for (double f : freq) {
            if (Double.isNaN(f)) {
                failure = true;
            }
        }
        
        // TODO: [ReconEP]:: the below shoudl be split into multiple helpers?
        if (!failure) {
        	boolean using_hc_clustering =false;
        	if (this.num_loci>40) {
        		using_hc_clustering =false;
        	}
        	
        	if (!using_hc_clustering) {
        	
	            boolean[] list_remove_haps = new boolean[freq.length];
	            int num_remove_hap = 0;
	            HashSet <String> aem_haps= new HashSet <String>();
	            
	            for (int h = 0; h < freq.length; h++) {
	                if (freq[h] < aem_parameters.aem_regional_cross_pool_freq_cutoff) {
	                    list_remove_haps[h] = true;
	                    num_remove_hap++;
	                } else {
	                	String tmp ="";
	                	for (int i=0;i< this.initial_Haps.global_haps_string[h].length;i++) {
	                		tmp=tmp+ this.initial_Haps.global_haps_string[h][i];
	                	}
	                	if (aem_haps.contains(tmp)) {
	                		list_remove_haps[h] = true;
	                        num_remove_hap++;
	                	}else {
	                		aem_haps.add(tmp);
	                		this.initial_Haps.global_haps_freq[h] = freq[h];
	                	}
	                    
	                }
	            }
	            
//	            if ( this.num_loci>35) {
//	            	aem_parameters.aem_hapset_size_min=2000;
//	            }
	            
	            double actual_cutoff = aem_parameters.aem_regional_cross_pool_freq_cutoff;
	            // If too many of them are below the regional frequency minimum,
	            // i.e., too few haps are remained 
	            if (aem_parameters.aem_hapset_size_min > freq.length - num_remove_hap) {
	                if(freq.length >= aem_parameters.aem_hapset_size_min) {
	                    list_remove_haps = Algebra.permute_sort_and_remove(freq.clone(), 
	                        aem_parameters.aem_hapset_size_min);
	                }
	            }
	            // Or, if too many haps are remained...
	            
	            
	            else if (aem_parameters.aem_hapset_size_max < freq.length - num_remove_hap) {
	                list_remove_haps = Algebra.permute_sort_and_remove(freq.clone(), 
	                    aem_parameters.aem_hapset_size_max);
	            }
	            
	            num_remove_hap = 0;
	            actual_cutoff=1;
	            for(int h=0; h < freq.length; h++) {
	                if(list_remove_haps[h]) {
	                    num_remove_hap++;
	                } else {
	                    this.initial_Haps.global_haps_freq[h] = freq[h];
	                    if(actual_cutoff > freq[h]) actual_cutoff = freq[h];
	                }
	            }
	            
	            
	            this.initial_Haps.remHaps(list_remove_haps, num_remove_hap);
	            System.out.println("Of the "
	                + freq.length
	                + " regional haplotypes, "
	                + num_remove_hap
	                + " were removed. The frequency cutoff was "
	                + actual_cutoff
	                + ".");           
	            this.final_Haps = this.initial_Haps;
	        }else {
	        	HierarchicalClustering hc= new HierarchicalClustering(  
	        				this.initial_Haps.global_haps_string, 
	        				freq, this.hc_cufoff, aem_parameters.aem_hapset_size_max, 
	        				this.initial_Haps);
	        		
	        	this.final_Haps = hc.final_hapconfig;
	        }
            
            
        }

    }
}
