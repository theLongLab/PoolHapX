package MiscFunctions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Properties;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import shapeless.newtype;
import spire.optional.intervalGeometricPartialOrder;


/*
 * The properties_file looks like this:
 * 
# PoolSimulator_SLiM Parameters
##########
# General: Commands and file locations
Input_Dir = D:\\PhD-Studying\\Informatics\\Project\\HIV_project\\PoolHapX_testing\\SLiM\\input
##fasta_Dir = /export/home/jhe/project/Viral_reconstruction/PoolHapX_test/PHX_perfect_20/input/fasta
Intermediate_Dir = D:\\PhD-Studying\\Informatics\\Project\\HIV_project\\PoolHapX_testing\\SLiM\\intermediate
Gold-Standard_Dir = D:\\PhD-Studying\\Informatics\\Project\\HIV_project\\PoolHapX_testing\\SLiM\\gold_standard
Slim_Output_Path = D:\\PhD-Studying\\Informatics\\Project\\HIV_project\\PoolHapX_testing\\SLiM\\input
Slim_Model = island_haploid
Proj_Name = 0_1
Is_Single_Population = false
Is_Ms_Output = false
DWGSIM = /gpfs/home/lmak/programs/DWGSIM-master/dwgsim
##########
Num_Pools = 9
Ref_Seq_Len = 9719
##########
Reference_Seq = D:\\PhD-Studying\\Informatics\\Project\\HIV_project\\PoolHapX_testing\\SLiM\\input\\HIV_HXB2.fa
Is_Perfect = false
Error_Rate_Per_Base = 0.0 
Coverage = 100
Read_Len = 100
Outer_Dist = 400
##########   
 */

public class PoolSimulator_SLiM {
	
    String input_dir;
    String inter_dir;
    String gs_dir;
    String vef_dir;
    String fastq_folder;
    String fasta_folder;
    String vef_folder;
    String slim_model;
    String slim_output_path;
    String project_name;
    Boolean is_single_population;
    // executables
    String msCMDLine ; 
    String slimCMDLine ;
    String dwgsimCMDLine ;
    
    // configurations set by the users
    int all_pool_haps;  // total number of haps in all pools.
    int haps_per_pool ;
    int[] haps_per_pop_arr;
    int num_pools ;
    int est_ind_pool ;
    double mutation_rate ;
    int num_var_pos ;
    int ref_seq_len ;
    String ref_seq_file_path; // full file path.     
    boolean is_perfect;
    double error_rate ;
    int coverage ;
    int read_len ;
    int outer_dist ;
    double hap_freq_cutoff;
    int weak_length;
    // intermediate global variables that are needed in the simulation.
    int actual_num_haps = 0; 
    int actual_num_vars = 0; 
    int[] sim_var_pos; 
    int[] hap2cts; // NOTE: This is also global count!
    double[] hap2allfreqs;  // # haps
    double[][] hap2infreqs; // # haps x # pools
    int[][] hap2varcomp;    // # haps x # vars 
    ArrayList<ArrayList<Integer>> hap2varpos= new ArrayList<ArrayList<Integer>>(); 
    // hap_id -> [alternate allele variant pos]
    HashMap<Integer, ArrayList<Integer>> pool2hapcomp = new HashMap<Integer, ArrayList<Integer>>();
    ArrayList<HashMap<String, Integer>> pool2allhapList=new ArrayList<HashMap<String, Integer>>();
    ArrayList<String> actual_hap_list = new ArrayList<String>();
    double var_burden_avg;
    boolean is_bacteria =false;
    int num_haps_pool;
    int random_var=1  ;
    
	public PoolSimulator_SLiM(String parameter_file) throws IOException {
	    InputStream is = new FileInputStream(parameter_file);
        Properties prop = new Properties();
        prop.load(is);
        
        this.input_dir = prop.getProperty("Input_Dir")+"/";
        this.inter_dir = prop.getProperty("Intermediate_Dir")+"/";
        this.gs_dir = prop.getProperty("Gold-Standard_Dir")+"/";
        new File(input_dir+"fasta/").mkdir();
        new File(input_dir+"fastq/").mkdir();
        new File(inter_dir+"vef/").mkdir();
        this.fasta_folder=this.input_dir + "fasta/";
        this.fastq_folder=this.input_dir + "fastq/";
        this.vef_folder=this.inter_dir+"vef/";
        this.slim_output_path=prop.getProperty("Slim_Output_Path")+"/";
        this.slim_model=prop.getProperty("Slim_Model");
        this.project_name=prop.getProperty("Proj_Name");
        this.is_single_population=Boolean.parseBoolean(prop.getProperty("Is_Single_Population"));
        this.dwgsimCMDLine = prop.getProperty("DWGSIM"); 
        this.num_pools = Integer.parseInt(prop.getProperty("Num_Pools"));
        this.ref_seq_len = Integer.parseInt(prop.getProperty("Ref_Seq_Len"));
        this.ref_seq_file_path = prop.getProperty("Reference_Seq"); 
        this.is_perfect = Boolean.parseBoolean(prop.getProperty("Is_Perfect"));
        if (is_perfect == false) {
            new File(input_dir+"sam/").mkdir();
            new File(input_dir+"bam/").mkdir();
            new File(input_dir+"vcf/").mkdir();
        }
        this.error_rate = Double.parseDouble(prop.getProperty("Error_Rate_Per_Base"));
        this.hap_freq_cutoff = Double.parseDouble(prop.getProperty("Hap_Freq_Cutoff"));
        this.coverage = Integer.parseInt(prop.getProperty("Coverage"));
        this.read_len = Integer.parseInt(prop.getProperty("Read_Len"));
        this.outer_dist = Integer.parseInt(prop.getProperty("Outer_Dist"));
        this.weak_length = Integer.parseInt(prop.getProperty("Weak_Length"));
        this.num_haps_pool = Integer.parseInt(prop.getProperty("Num_Haps_Pool"));
        is.close();
 
	}
	
	/**
	 * // Step 1A: Processing slim standard outcome
	 * Figure out:
	 * i) the number of types of haplotypes and 
	 * ii) the non-degenerate variant positions.
	 * 
	 * @throws IOException
	 */
	public void processing_standard_outcome(Boolean is_single_population) throws IOException{
		System.out.println("Step 1A: Figure out i) the number of types of haplotypes and ii) "
				    + "the non-degenerate variant positions.\n");
		BufferedReader br = new BufferedReader(new FileReader(gs_dir + project_name + "_"+slim_model+".out")); 
		ArrayList<String> index2varpos = new ArrayList<String>();
		String currLine = br.readLine(); // header
		String[] tmpcurrpos = currLine.split(" "); 
		if(is_single_population==true) {
			HashMap<String, Integer> hapforpool= new HashMap<String, Integer>();
			this.pool2allhapList.add(hapforpool);
			while(!tmpcurrpos[0].equals("Mutations:")) { // Read those lines before "Mutations:"
				currLine = br.readLine();
				tmpcurrpos = currLine.split(" ");
			}
			currLine = br.readLine();
			tmpcurrpos = currLine.split(" ");
		}else {
			while(!tmpcurrpos[0].equals("Populations:")) { // Read those lines before "Populations:"
				currLine = br.readLine();
				tmpcurrpos = currLine.split(" ");
			}
			currLine = br.readLine();
			while(!tmpcurrpos[0].equals("Mutations:")) { // Read those lines before "Mutations:"
				HashMap<String, Integer> hapforpool= new HashMap<String, Integer>();
				this.pool2allhapList.add(hapforpool);
				currLine = br.readLine();
				tmpcurrpos = currLine.split(" ");
			}
			currLine = br.readLine();
			tmpcurrpos = currLine.split(" ");
		}
		// get rid of variant positions that occur at beginning of weak_length
		// and at end of the weak_length
		//int weak_length = this.read_len*2+this.outer_dist;
		int front_cutoff = weak_length;
		int end_cutoff = this.ref_seq_len -weak_length;
		ArrayList<String> var_pos_list = new ArrayList<String>(); 
		while(!tmpcurrpos[0].equals("Genomes:")&&!tmpcurrpos[0].equals("Individuals:")) {
			//Read those lines under the "Mutations" section
			if(Integer.parseInt(tmpcurrpos[3])>front_cutoff && Integer.parseInt(tmpcurrpos[3])<end_cutoff ) {
				if(!var_pos_list.contains(tmpcurrpos[3])) {
					var_pos_list.add(tmpcurrpos[3]);
					index2varpos.add(tmpcurrpos[0]+"_"+tmpcurrpos[3]);
				}
			}	
			currLine = br.readLine();
			tmpcurrpos = currLine.split(" ");
		}
		System.out.println("Index to variants positions: "+ index2varpos);
		this.num_var_pos=index2varpos.size();
		this.sim_var_pos=new int[num_var_pos];
		int[] sim_var_index=new int[num_var_pos];
		for(int p=0;p<num_var_pos;p++) {
			String[] index2varposarray = index2varpos.get(p).split("_");
			sim_var_index[p]=Integer.parseInt(index2varposarray[0]);
			sim_var_pos[p]=Integer.parseInt(index2varposarray[1]);

		}
		// sort both sim_var_index[] and sim_var_pos[] at the same time according to the sim_var_pos[]
		for(int i = 0; i < sim_var_pos.length - 1;i++) {
			for(int j = 0; j < sim_var_pos.length - 1 - i; j++) {
				if(sim_var_pos[j] >sim_var_pos[j+1]) {
					int temp_pos = sim_var_pos[j];
					sim_var_pos[j] = sim_var_pos[j+1];
					sim_var_pos[j+1] = temp_pos;
					int temp_index = sim_var_index[j];
					sim_var_index[j] = sim_var_index[j+1];
					sim_var_index[j+1] = temp_index;
				}
			}
		}

		if(!tmpcurrpos[0].equals("Genomes:")) {
			while(!tmpcurrpos[0].equals("Genomes:")) {
				 currLine = br.readLine();
				 tmpcurrpos = currLine.split(" ");
			}
		}
		currLine = br.readLine();
		// pool2allhapList: haplotypes and their corresponding number within each pool
		// hapsHS:all the haplotypes and their corresponding number cross all pools
		// For single_population, the pool2allhapList and hapsHS will be the same,
		// we just need to use the hapsHS and randonly assign those haplotypes to different pools
		// For island model, we generate the gold_standard and fastq file using pool2allhapList
		HashMap<String, Integer> hapsHS = new HashMap<String, Integer>();
		int num_hap=0;
		int max_num_haps = this.num_pools* this.num_haps_pool;
		while(currLine!=null) {
			tmpcurrpos = currLine.split(" ");
			int curr_pool_index= Integer.parseInt(tmpcurrpos[0].split(":")[0].split("p")[1])-1;
			if ((curr_pool_index<= this.num_pools) || (is_single_population) ){
				if(tmpcurrpos.length>2) {
					num_hap++;
					
					ArrayList<Integer> currvarpos = new ArrayList<Integer>();
					String curr_hap="";
					for (int p=2;p<tmpcurrpos.length;p++){
						currvarpos.add(Integer.parseInt(tmpcurrpos[p]));
					}
					for(int p=0;p<num_var_pos;p++) {
						int curr_var_index=sim_var_index[p];
						if(currvarpos.contains(curr_var_index)) {
							curr_hap=curr_hap+1;
						}else {
							curr_hap=curr_hap+0;
						}
					}
					if (!hapsHS.containsKey(curr_hap)) {
						hapsHS.put(curr_hap, 1);
					}else {
						int tmpCt =  hapsHS.get(curr_hap) + 1;
						hapsHS.put(curr_hap, tmpCt);
					}
	
					if(!pool2allhapList.get(curr_pool_index).containsKey(curr_hap)) {
						this.pool2allhapList.get(curr_pool_index).put(curr_hap, 1);
					}else if(pool2allhapList.get(curr_pool_index).containsKey(curr_hap)) {
						int tmpCt = pool2allhapList.get(curr_pool_index).get(curr_hap)+1;
						this.pool2allhapList.get(curr_pool_index).put(curr_hap, tmpCt);
					}
				}
			}
			currLine = br.readLine();
		}
		br.close();
		
		
		
		this.all_pool_haps=num_hap;
		double num_hap_cutoff = (double)all_pool_haps*this.hap_freq_cutoff; 
		ArrayList<String> low_freq_hap = new ArrayList<String>();
		for (String h : hapsHS.keySet()) {
			if(hapsHS.get(h)<num_hap_cutoff) {
				low_freq_hap.add(h);
				
			}
		}
		System.out.println(low_freq_hap);
		for(int i=0;i<low_freq_hap.size();i++) {
			all_pool_haps=all_pool_haps-hapsHS.get(low_freq_hap.get(i));
			hapsHS.remove(low_freq_hap.get(i));
			for(int pool=0;pool<pool2allhapList.size();pool++) {
				if(pool2allhapList.get(pool).containsKey(low_freq_hap.get(i))) {
				pool2allhapList.get(pool).remove(low_freq_hap.get(i));
				}
			}
		}
		
		this.actual_num_haps = hapsHS.size();
        this.hap2varcomp = new int[actual_num_haps][num_var_pos]; 
        this.hap2cts = new int[actual_num_haps]; 
        int hap = 0; 
        double var_burden_ct = 0.0; 
        int[] true_var_pos = new int[num_var_pos];
        for (String h : hapsHS.keySet()) {
        	this.actual_hap_list.add(h);
            String[] tmpHapComp = h.split("");
            this.hap2varpos.add(new ArrayList<Integer>());
            for (int p = 0; p < num_var_pos; p++) {
                int tmpAllele = Integer.parseInt(tmpHapComp[p]); 
                this.hap2varcomp[hap][p] = tmpAllele; 
                if (tmpAllele == 1) {
                    true_var_pos[p] = 1;    
                    // If this variant position is represented by at least one alternate allele, 
                    // then it's a true variant position.
                    this.hap2varpos.get(hap).add(sim_var_pos[p]);
                }
                var_burden_ct += (double) tmpAllele; 
            }
            hap2cts[hap] = hapsHS.get(h);
            hap++; 
        }
		System.out.println("Haplotypes per pools list: " + pool2allhapList);
        this.actual_num_vars = SimpleMath.sum(true_var_pos); 
        this.var_burden_avg = var_burden_ct / (double) actual_num_haps; 
	}
	
	/**
	 * Step 1A: Processing slim msoutcome
	 */
	public void processing_ms_outcome() throws IOException{
	  System.out.println("Step 1A: Figure out i) the number of types of haplotypes and ii) "
		    + "the non-degenerate variant positions.\n");
	    BufferedReader br = new BufferedReader(new FileReader(gs_dir + project_name + "_slim.txt")); 
	    String currLine = br.readLine(); // header
	    String[] tmpVarPos = currLine.split(" ");
	    while(!tmpVarPos[0].equals("#OUT:")) {
	    	currLine = br.readLine();
	    	tmpVarPos = currLine.split(" ");
	    }
	    this.all_pool_haps = Integer.parseInt(tmpVarPos[tmpVarPos.length-1]);
	    this.haps_per_pool = all_pool_haps/this.num_pools;
	    
	    while(!tmpVarPos[0].equals("segsites:")) {
	    	currLine = br.readLine();
	    	tmpVarPos = currLine.split(" ");
	    }
	    this.num_var_pos = Integer.parseInt(tmpVarPos[1]);
	    this.sim_var_pos = new int[num_var_pos];

	    currLine = br.readLine();	    
	    tmpVarPos = currLine.split(" "); // positions
	    for (int p = 1; p <= num_var_pos; p++) {
	        sim_var_pos[p - 1] = (int) Math.floor(Double.parseDouble(tmpVarPos[p])
	         	* ref_seq_len);
	         // If there isn't enough of a difference between adjacent fractions generated by ms.
	        if (p > 1 && sim_var_pos[p - 1] == sim_var_pos[p - 2]) {
	             sim_var_pos[p - 1]++;
	        }
	    }
        currLine = br.readLine();
        HashMap<String, Integer> hapsHS = new HashMap<String, Integer>(); 
        for (int h = 0; h < all_pool_haps; h++) {
            if (!hapsHS.containsKey(currLine)) hapsHS.put(currLine, 0);
            int tmpCt =  hapsHS.get(currLine) + 1;
            hapsHS.put(currLine, tmpCt);
            currLine = br.readLine();
        }
        br.close();
        this.actual_num_haps = hapsHS.size();
        this.hap2varcomp = new int[actual_num_haps][num_var_pos]; 
        this.hap2cts = new int[actual_num_haps]; 
        int hap = 0; 
        double var_burden_ct = 0.0; 
        int[] true_var_pos = new int[num_var_pos];
        for (String h : hapsHS.keySet()) {
            String[] tmpHapComp = h.split("");
            this.hap2varpos.add(new ArrayList<Integer>());
            for (int p = 0; p < num_var_pos; p++) {
                int tmpAllele = Integer.parseInt(tmpHapComp[p]); 
                this.hap2varcomp[hap][p] = tmpAllele; 
                if (tmpAllele == 1) {
                    true_var_pos[p] = 1;    
                    // If this variant position is represented by at least one alternate allele, 
                    // then it's a true variant position.
                    this.hap2varpos.get(hap).add(sim_var_pos[p]);
                }
                var_burden_ct += (double) tmpAllele; 
            }
            hap2cts[hap] = hapsHS.get(h);
            hap++; 
        }
        this.actual_num_vars = SimpleMath.sum(true_var_pos); 
        this.var_burden_avg = var_burden_ct / (double) actual_num_haps; 
	}

	
	/**
	 * Step 1B: Report properties of the simulated haplotypes to the user 
	 * to check if they're acceptable.
	 * @param prefix
	 */
	public void ms_reports() throws IOException {
        System.out.println("Step 1B: Report properties of the simulated haplotypes to the user "
            + "to check if they're acceptable.");
        
        int[] pwDifference = new int[actual_num_haps * (actual_num_haps - 1) / 2];
        int compare = 0; 
        this.hap2allfreqs = new double[actual_num_haps]; 
        for (int h = 0; h < actual_num_haps; h++) {
            for (int i = h + 1; i < actual_num_haps; i++) {
                for (int p = 0; p < num_var_pos; p++)
                    if (hap2varcomp[h][p] != hap2varcomp[i][p]) 
                        pwDifference[compare]++;
                if (pwDifference[compare] == 0) System.out.println(h + "\t" + i + "\t");
                // System.out.print(pwDifference[compare] + "\t");
                compare++; 
            }
            hap2allfreqs[h] = (double) hap2cts[h] / all_pool_haps;
        }
        Arrays.sort(pwDifference);
        double meanPWDiff = SimpleMath.mean(pwDifference);
        double stdPWDiff = SimpleMath.stdev(pwDifference);
        int[] sortedCts = new int[hap2cts.length];
        System.arraycopy(hap2cts, 0, sortedCts, 0, hap2cts.length);
        Arrays.sort(sortedCts);
        double meanCts = SimpleMath.mean(hap2cts);
        double stdCts = SimpleMath.stdev(hap2cts);
        
        System.out.println("There are " + actual_num_haps + " across-pool haplotypes.");
        System.out.println("There are " + actual_num_vars + " true variant positions.");
        System.out.println("There average mutational burden per haplotype is " 
        + var_burden_avg + " variants.");
        System.out.println("The average pairwise difference is " + meanPWDiff 
            + " and the standard deviation is " + stdPWDiff + ".");
        System.out.println("The minimum pairwise difference is " + pwDifference[0] 
            + " and the maximum is " + pwDifference[pwDifference.length - 1] + "."); 
        System.out.println("The average all-pool count per haplotype is " 
            + meanCts + " and the standard deviation is " + stdCts + ".");
        System.out.println("The minimum all-pool count is " + sortedCts[0] 
            + " and the maximum is " + sortedCts[sortedCts.length - 1] + "."); 
        
        PrintWriter pw = new PrintWriter(new FileWriter(gs_dir 
            + "PD.simulation_summary.txt", true));   // gs_dir/c.simulation_summary.txt
        pw.append("Project_Name"+"\t"+"Total_Hap_Count"+"\t"+"Num_Var_Pos"+"\t"
            +"Ave_Mutation_Burden_Per_Hap"+"\t"+"Ave_Pairwise_Diff"+"\t"+"Std"
        		+"\t"+"Min_Pairwise_Diff"+"\t"+"Max_Pairwise_Diff"+"\t"
            +"Ave_Count_Per_Hap"+"\t"+"Std_Per_Hap"+"\t"+"Min_Count_Per_Hap"
        		+"\t"+"Max_Count_Per_Hap"+"\n");
        pw.append(project_name + "\t" + actual_num_haps + "\t" + actual_num_vars 
            + "\t" + var_burden_avg + "\t" + meanPWDiff + "\t" + stdPWDiff + "\t" 
            + pwDifference[0] + "\t" + pwDifference[pwDifference.length - 1] + "\t" + meanCts  
            + "\t" + stdCts + "\t" + sortedCts[0] + "\t" + sortedCts[sortedCts.length - 1] + "\n"+"\n");
        pw.close();
        // System.out.print("Is this acceptable? ");
        // answer = reader.next();
        // reader.close();
    // } while (!answer.equals("Y"));   // Basically, simulate haplotypes until the distribution makes me happy.
	}
            
    /**
     * Step 2: Assign each haplotype individual to a patient, 
     * and write all of the gold standard files for PoolHapX.
     * @param prefix
     * @throws IOException
     * @throws InterruptedException
     */
	public void assign_haps_to_pools(Boolean is_single_population) throws IOException {                          
        System.out.println("\nStep 2: Assign each haplotype individual to a patient, "
            + "and write all of the gold standard files.\n");
        this.hap2infreqs = new double[actual_num_haps][num_pools];
        int[][] var2incts = new int[actual_num_vars][num_pools];
    	double[][] var2infreqs = new double[actual_num_vars][num_pools];
        if(is_single_population==true) {
        	int[][] hap2incts = new int[actual_num_haps][num_pools];
        	int random_count =0;
            for (int h = 0; h < actual_num_haps; h++) {
                while (hap2cts[h] != 0) {
                	Random seedRandom = new Random( random_count*this.random_var);
                	this.random_var++;
                	random_count++;
                    int currPool = seedRandom.nextInt( num_pools);
//                    System.out.println(h+"\t"+currPool);
                    if (!this.pool2hapcomp.containsKey(currPool)) {
                    	this.pool2hapcomp.put(currPool, new ArrayList<Integer>());
                    	this.pool2hapcomp.get(currPool).add(h);
                    	hap2incts[h][currPool]++; 
                    	for (int v = 0; v < num_var_pos; v++) var2incts[v][currPool] += hap2varcomp[h][v];
                        hap2cts[h]--;
                    }else {
                    	this.pool2hapcomp.get(currPool).add(h);
                    	hap2incts[h][currPool]++; 
                        for (int v = 0; v < num_var_pos; v++) var2incts[v][currPool] += hap2varcomp[h][v];
                        hap2cts[h]--; 
                    }
                }
            }
            for(int h=0;h<actual_num_haps;h++) {
            	for(int p = 0; p < num_pools; p++) {
            		
            		int haps_this_pool = pool2hapcomp.get(p).size();
                	hap2infreqs[h][p] = (double) hap2incts[h][p] /haps_this_pool;
            	}
            }
            for(int p = 0; p < num_pools; p++) {
            	int haps_this_pool = pool2hapcomp.get(p).size();
                for(int v = 0; v < actual_num_vars; v++) {
                    var2infreqs[v][p] = (double) var2incts[v][p] / haps_this_pool;
                }
            }
            
            
        }else if(is_single_population==false) {
        	if(this.num_pools!=pool2allhapList.size()) {
        		System.out.println("Number of pools in property file is not consistent with the SLiM output");
        	}
        	this.haps_per_pop_arr= new int[this.num_pools];
        	for(int p=0;p<haps_per_pop_arr.length;p++) {
        		int curr_num_hap=0;
        		for(String h:pool2allhapList.get(p).keySet()) {
        			curr_num_hap=curr_num_hap+pool2allhapList.get(p).get(h);
        		}
        		haps_per_pop_arr[p]=curr_num_hap;
        	}
        	// Generate hap2infreqs according to the pool2allhapList
        	for(int h=0;h<actual_num_haps;h++) {
        		String curr_hap = actual_hap_list.get(h);
        		for(int p=0;p<num_pools;p++) {
        			if(pool2allhapList.get(p).containsKey(curr_hap)) {
        				double hap2poolfre = (double)pool2allhapList.get(p).get(curr_hap)/
        						(double)haps_per_pop_arr[p];
        				hap2infreqs[h][p]=hap2poolfre;
        			}else {
        				hap2infreqs[h][p]=0.0;
        			}
        		}
        	}
        	// Generate var2infreqs according to the pool2allhapList
        	for(int p=0;p<num_pools;p++) {
        		for(String h:pool2allhapList.get(p).keySet()) {
        			for(int v=0;v<sim_var_pos.length;v++) {
        				if(h.split("")[v].equals("1")) {
        					var2incts[v][p]=var2incts[v][p]+pool2allhapList.get(p).get(h);
        				}
        			}
        		}
        	}
        	for(int p=0;p<num_pools;p++) {
        		for(int v=0;v<sim_var_pos.length;v++) {
        			var2infreqs[v][p] = (double)var2incts[v][p]/(double)haps_per_pop_arr[p];
        		}
        	}
        } // end of if
        
        	BufferedWriter bw_inter = new BufferedWriter(new FileWriter(gs_dir + project_name + 
                "_haps.inter_freq_vars.txt"));
        	bw_inter.write("Hap_ID");
            for(int h = 0; h < actual_num_haps; h++)
            	bw_inter.write("\t" + "h"+h);
            bw_inter.write("\nFreq");
            for(int h = 0; h < actual_num_haps; h++)
            	bw_inter.write("\t" + hap2allfreqs[h]);
            bw_inter.write("\n");
            for(int v = 0; v < actual_num_vars; v++){
            	bw_inter.write("0;" + sim_var_pos[v] + ";" + sim_var_pos[v] + ";0:1");
                for(int h = 0; h < actual_num_haps; h++)
                	bw_inter.write("\t" + hap2varcomp[h][v]);
                bw_inter.write("\n");
            } bw_inter.close();
            
            
//            
//            for (int p = 0; p < num_pools; p++) {
//            	int haps_this_pool = pool2hapcomp.get(p).size();
//                PrintWriter pw = new PrintWriter(fasta_folder + project_name + "_p" + p + ".fa");
//                	for (int h = 0; h < haps_this_pool; h++) {
//                		int currHap = this.pool2hapcomp.get(p).get(h);
//                		pw.append(">Haplotype_" + currHap + " \n");
//                		for (String s : allSimHaps[currHap]) pw.append(s);
//                		pw.append("\n\n");
//                	}
//                pw.close();
//            	}
            
            if (is_bacteria) {	
            	
//            	fastq_folder + project_name + "_p" + p, 
            	
            	for (int p = 0; p < num_pools; p++) {
            		
            		int haps_this_pool = pool2hapcomp.get(p).size();
            		int[] hap_count = new int[all_pool_haps];
            		for (int i=0;i< hap_count.length;i++) {
            			hap_count[i] =0;
            		}
            		for(int h = 0; h < haps_this_pool; h++) {
            			hap_count[this.pool2hapcomp.get(p).get(h)]++;
            		}
		        	BufferedWriter bw_bed = new BufferedWriter(new FileWriter(fasta_folder + project_name + "_p" + 
		        			Integer.toString(p) + ".bed"));
		        	
		        	for(int h = 0; h < haps_this_pool; h++) {
		        		for(int t = 0; t < hap_count[h]; t++) {
		        			for(int v = 0; v < actual_num_vars; v++){
		        			
		        				String Chrom = "Haplotype_" + Integer.toString(h)+"_"+Integer.toString(t);
		        				bw_bed.write(Chrom+"\t" +  (sim_var_pos[v]-750)+ "\t"+ 
				        				(sim_var_pos[v]+750)+"\n");
		        			}
		        		}
		        	}
		        		

			        	
			       bw_bed.close();
            	}
            }
        	
  
            

            BufferedWriter bw = new BufferedWriter(new FileWriter(gs_dir + project_name + "_haps.intra_freq.txt"));
            bw.write("Hap_ID");
            for(int h = 0; h < actual_num_haps; h++)
                bw.write("\t" + "h"+h);
            bw.write("\n");
            for(int p = 0; p < num_pools; p++){
                bw.write(project_name + "_p"+p);
                for(int h = 0; h < actual_num_haps; h++)
                    bw.write("\t" + hap2infreqs[h][p]);
                bw.write("\n");
            }
            bw.close();

            if(is_perfect == true) {
            	bw = new BufferedWriter(new FileWriter(inter_dir + project_name + "_vars.intra_freq.txt"));
            }else {
            	bw = new BufferedWriter(new FileWriter(gs_dir + project_name + "_vars.intra_freq.txt"));
            }
            bw.write("Pool_ID");
            for (int p = 0; p < num_pools; p++)
                bw.write("\t" + project_name + "_p" + p); 
            bw.write("\n");
            for (int v = 0; v < actual_num_vars; v++){
                bw.write("0;" + sim_var_pos[v] + ";" + sim_var_pos[v] + ";0:1");
                for (int p = 0; p < num_pools; p++)
                    bw.write("\t" + var2infreqs[v][p]);
                bw.write("\n");
            } bw.close();
	}
	

	/**
	 * Step 3: Make all of the patient FASTA files, and simulate reads for them.Generate FASTQ files  
     * 
	 * @throws IOException
	 * @throws InterruptedException
	 */
        // 
	public void generate_fastq(Boolean is_single_population)  throws IOException, InterruptedException {
        System.out.println("Step 3: Make all of the patient FASTA files.");
        BufferedReader br = new BufferedReader(new FileReader(ref_seq_file_path));
        String[] refSequence = new String[ref_seq_len];
        String currLine = br.readLine();
        int i = 0;
        while (currLine != null) {
            if (currLine.contains(">")) {
                currLine = br.readLine();
                continue;
            }            
            for (String b : currLine.split("")) {
                refSequence[i] = b;
                i++;
            }            
            currLine = br.readLine();
        }
        br.close();
        // Step 5a) Simulate single nucleotide polymorphisms on the reference sequence.
        System.out.println(
            "\nStep 3A Simulate single nucleotide polymorphisms on the reference sequence.\n");
        BufferedWriter bw = new BufferedWriter(new FileWriter(gs_dir + project_name + "_mutations.txt"));
        bw.append("Variants"+"\t"+"Positions"+"\t"+"Ref_Allele"+"\t"+"Alt_Allele"+"\n");
        String[] allAltAlleles = new String[actual_num_vars];
        for (int v = 0; v < actual_num_vars; v++) {
            allAltAlleles[v] = simVariant(refSequence[sim_var_pos[v] - 1]);
            bw.append(v + "\t" + sim_var_pos[v] + "\t" + refSequence[sim_var_pos[v] - 1] 
                + "\t" + allAltAlleles[v] + "\n");
        }
        bw.close();        
        
     // Step 5b) Add simulated mutations to finish the full-length simulated haplotypes.
        System.out.println(
            "\nStep 3B Add simulated mutations to finish the full-length simulated haplotypes.\n");

        String[][] allSimHaps = new String[actual_num_haps][ref_seq_len];
        for (int h = 0; h < actual_num_haps; h++) {
            for (int p = 0; p < ref_seq_len; p++) {
                int pos = find(sim_var_pos, p + 1);
                if (pos != -1) {
                    if (hap2varcomp[h][pos]==1) {
                        allSimHaps[h][p] = allAltAlleles[pos];
                    } else {
                        allSimHaps[h][p] = refSequence[p];
                    }                    
                } else {
                    allSimHaps[h][p] = refSequence[p];
                }
            }
        }
        
        // Simulate all of the pool FastA and FastQ files, given the distribution of haplotypes in step 3.
        // HashMap<Integer, ArrayList<Integer>> pool2hapcomp = new HashMap<Integer, ArrayList<Integer>>(); // pool id -> [hap_ids]
        System.out.println("\nStep 3C: Simulate all of the pool FastA files, given the distribution"
            + " of haplotypes in step 3.\n");
        
        if(is_single_population==true) {
        	for (int p = 0; p < num_pools; p++) {
        	int haps_this_pool = pool2hapcomp.get(p).size();
        	HashMap<Integer, Integer> hap_count = new HashMap<Integer, Integer>();
            PrintWriter pw = new PrintWriter(fasta_folder + project_name + "_p" + p + ".fa");
            	for (int h = 0; h < haps_this_pool; h++) {
            		int currHap = this.pool2hapcomp.get(p).get(h);
            		if (!hap_count.containsKey(currHap)) {
            			hap_count.put(currHap, 0);
	            		pw.append(">Haplotype_" + currHap + "_0 \n");
	            		for (String s : allSimHaps[currHap]) pw.append(s);
	            		pw.append("\n\n");
            		}else {
            			hap_count.put(currHap , hap_count.get(currHap)+1);
            			pw.append(">Haplotype_" + currHap +  "_"+ hap_count.get(currHap)+" \n");
	            		for (String s : allSimHaps[currHap]) pw.append(s);
	            		pw.append("\n\n");
            		}
            	}
            pw.close();
        	}
        }else if(is_single_population==false) {
        	for (int p = 0; p < num_pools; p++) {
                PrintWriter pw = new PrintWriter(fasta_folder + project_name + "_p" + p + ".fa");
                for(int h=0;h < actual_num_haps;h++) {
                	String curr_hap = actual_hap_list.get(h);
                	if(pool2allhapList.get(p).containsKey(curr_hap)) {
                		for(int currhap_num=0;currhap_num<pool2allhapList.get(p).get(curr_hap);currhap_num++) {
                    		pw.append(">Haplotype_" + h   + " \n");
                    		for (String s : allSimHaps[h]) pw.append(s);
                    		pw.append("\n\n");
                		}
                	}
                }
                pw.close();
        	}
        }
        
        System.out.println("\nStep 3D: Simulate all of the pool FastQ files, given the distribution"
                + " of haplotypes in step 3.\n");
        
        for (int p = 0; p < num_pools; p++) {
        	int haps_this_pool;
        	if(is_single_population==true) {
        		haps_this_pool = pool2hapcomp.get(p).size();
        	}else {
        		haps_this_pool=haps_per_pop_arr[p];
        	}
        	System.out.println("Haplotypes within this pool:" + haps_this_pool);
            double coverage_set = (double)coverage/(double)haps_this_pool;
            
           if (is_bacteria) {
	            ProcessBuilder CMDLine = new ProcessBuilder(dwgsimCMDLine,
	                fasta_folder + project_name + "_p" + p + ".fa", 
	                fastq_folder + project_name + "_p" + p, 
	                "-e", Double.toString(error_rate),
	                "-E", Double.toString(error_rate), "-C", Double.toString(coverage_set),
	                "-1", Integer.toString(read_len),
	                "-2", Integer.toString(read_len),
	                "-r", "0",
	                "-F", "0",
	                "-H",
	                "-d", Integer.toString(outer_dist),
	                "-o", "1",
	                "-s", "0",
	                "-y", "0",
	                "-x", fasta_folder + project_name + "_p" + 
		        			Integer.toString(p) + ".bed");
	            
	            Process CMDProcess = CMDLine.start();
	            CMDProcess.waitFor();
	            System.out.println("Finished simulating reads for pool " + p + ".");
           }else {
        	   System.out.println(dwgsimCMDLine+ fasta_folder + project_name + "_p" + p + ".fa"
        			   +   fastq_folder + project_name + "_p" + p+ " -e "+ Double.toString(error_rate)
        			   +  "-E "+ Double.toString(error_rate)+ " -C "+  Double.toString(coverage_set)
        			   +  "-1"+  Integer.toString(read_len)+ "-2"+  Integer.toString(read_len)+ 
        			   " -r "+  "0"+ " -F "+  "0"+  " -H "+ " -d "+  Integer.toString(outer_dist)+
        			   " -o "+  "1"+ " -s "+ "0"+ " -y "+  "0"
      	                );
        	   
        	   ProcessBuilder CMDLine = new ProcessBuilder(dwgsimCMDLine,
   	                fasta_folder + project_name + "_p" + p + ".fa", 
   	                fastq_folder + project_name + "_p" + p, 
   	                "-e", Double.toString(error_rate),
   	                "-E", Double.toString(error_rate), "-C", Double.toString(coverage_set),
   	                "-1", Integer.toString(read_len),
   	                "-2", Integer.toString(read_len),
   	                "-r", "0",
   	                "-F", "0",
   	                "-H",
   	                "-d", Integer.toString(outer_dist),
   	                "-o", "1",
   	                "-s", "0",
   	                "-y", "0");
   	            
   	            Process CMDProcess = CMDLine.start();
   	            CMDProcess.waitFor();
   	            System.out.println("Finished simulating reads for pool " + p + ".");
           }
        }
	}
	
	/**
	 * Step 4: Convert each patient's FASTQ file(s) to a VEF file. 
	 * @param 
	 * @return
	 */
	
	public void generate_VEF() throws IOException, InterruptedException {
        System.out.println("Concurrently, Step 4: Convert each patient's FASTQ file(s) to a VEF file. \n");
	    int startOne = 0, startTwo = 0, endOne = 0, endTwo = 0; 
	    for(int p=0;p<num_pools;p++) {
	        ProcessBuilder CMDLine = new ProcessBuilder("gunzip", fastq_folder + project_name 
	            + "_p" + p + ".bwa.read1.fastq.gz");
            Process CMDProcess2 = CMDLine.start();  
            CMDProcess2.waitFor();
            
            BufferedReader br = new BufferedReader(new FileReader(fastq_folder + project_name + 
                "_p" + p + ".bwa.read1.fastq"));
            String currLine = br.readLine();
            BufferedWriter bw  = new BufferedWriter(new FileWriter(vef_folder + project_name  
                + "_p" + p + ".vef")); 
            while (currLine != null) {
                String[] readInfo = currLine.split("_");
                StringBuilder readOutput = new StringBuilder(); 
                readOutput.append(currLine.trim() + "\t"); 
                int currID = Integer.parseInt(readInfo[1]); 
                int readOne = Integer.parseInt(readInfo[2]);
                int readTwo = Integer.parseInt(readInfo[3]); 
                if (readOne < readTwo)  {   // Organizes the reporting of the last, position-reporting columns of the VEF file properly.
                    startOne = readOne;
                    startTwo = readTwo; 
                } else {
                    startOne = readTwo;
                    startTwo = readOne;
                }
                endOne = startOne + read_len - 1; // Last included base in the read.
                endTwo = startTwo + read_len - 1; // Last included base in the read.
                Boolean varPosInRead = false;
                for (int v = 0; v < actual_num_vars; v++) {
                    if ((startOne <= sim_var_pos[v] && sim_var_pos[v] <= endOne) || 
                        (startTwo <= sim_var_pos[v] && sim_var_pos[v] <= endTwo)) {
                        if (this.hap2varpos.get(currID).contains(sim_var_pos[v])) {
                            readOutput.append(sim_var_pos[v] + "=1;");
                        } else {
                            readOutput.append(sim_var_pos[v] + "=0;");
                        } varPosInRead = true;
                    }
                }
                for (int l = 0; l < 4; l++) currLine = br.readLine(); 
                // Skip the next three lines (bases, +, base quality) and take the next name
                if (!varPosInRead) continue; 
                // If this read does not contain any segregating sites, it is not included in the final VEF..
                readOutput.append("\t//\t" + startOne + "\t" + endOne  + "\t" 
                    + startTwo + "\t" + endTwo  + "\n");  
                // If it does, then we need to output the variant info into the pool VEF.
                bw.write(readOutput.toString());
            }
            br.close();
            bw.close();
            System.out.println("Finished converting FASTQ format to VEF format for pool " + p + ".");
            
            // gzip the gastq file back.
            ProcessBuilder CMDLine_gzip = new ProcessBuilder("gzip", fastq_folder + project_name 
                + "_p" + p + ".bwa.read1.fastq");
            Process CMDProcess_gzip = CMDLine_gzip.start();  
            CMDProcess_gzip.waitFor();
	    }
	}
	
	static String simVariant(String refBase) {
	        ArrayList<String> bases = new ArrayList<String>();
	        if (!refBase.equals("A")) bases.add("A");
	        if (!refBase.equals("C")) bases.add("C");
	        if (!refBase.equals("G")) bases.add("G");
	        if (!refBase.equals("T")) bases.add("T");
	        Random seedRandom = new Random( 0);
	        return bases.get(seedRandom.nextInt(3));
	}

	static int find(int[] a, int target) {
	        int index = Arrays.binarySearch(a, target);
	        return (index < 0) ? -1 : index;
	}
	  
	public static void main(String[] args) throws IOException, InterruptedException {
		String parameter= args[0];//"D:\\PhD-Studying\\Informatics\\Project\\HIV_project\\PoolHapX_testing\\SLiM\\input\\PoolSimulator.properties";//
		InputStream is = new FileInputStream(parameter);
        Properties prop = new Properties();
        prop.load(is);
        Boolean is_perfect = Boolean.parseBoolean(prop.getProperty("Is_Perfect"));
        Boolean is_single_population = Boolean.parseBoolean(prop.getProperty("Is_Single_Population"));
        Boolean is_ms_output = Boolean.parseBoolean("Is_Ms_Output");
        is.close();
        //1st: Read the property file
		PoolSimulator_SLiM ps=new PoolSimulator_SLiM(parameter);
		//2nd:Processing SLiM output
		if(is_ms_output==false) {
			ps.processing_standard_outcome(is_single_population);
		}else {
			ps.processing_ms_outcome();
		}
		//3rd: Report the properties of the simulated haplotypes
		ps.ms_reports();
	    //4th: Assign haplotypes to individuals
	    ps.assign_haps_to_pools(is_single_population);
	    //5th step: Generate FASTA and FASTQ files
	    ps.generate_fastq(is_single_population);
	    //6th step: Only if is_perfect=true, run the 6th step
	    if(is_perfect==true) {
	    	ps.generate_VEF();
	    }
	}
}