package MiscFunctions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashSet;

public class GC1_CompareHaps {
	
	/*
	 * generate_inter_file()
	 * This function will generate one/two files under the output_dir.
	 * The "all" means to include all possible haplotypes and ignoring the "?".
	 * A file_name end up with "_all.inter_freq_haps.txt" will be generated.
	 * The "complete" means to include those haplotypes without "?" only.
	 * If the complete haplotypes do exits, a file end up with 
	 * "_complete.inter_freq_haps.txt" will be generated. Otherwise, it will
	 * print message " no complete_inter_file is generate".
	 */
	public static void generate_inter_file(String gcf_file,String vars_intra_file,
			String all_inter_file, String complete_inter_file) 
			throws IOException, InterruptedException  {
		
		ArrayList<String> loci_string_list=new ArrayList<String>();
		//hap_full_string_list ignoring "?", and take into account all possible haploypes
		ArrayList<String> hap_all_string_list=new ArrayList<String>();	
		ArrayList<ArrayList<String>> all_hap_seq_listlist=new ArrayList<ArrayList<String>>();
		//patial_hap_seq_listlist only takes those complete haplotypes into account(those without ?)
		ArrayList<ArrayList<String>> complete_hap_seq_listlist=new ArrayList<ArrayList<String>>();
		int complete_hap=0;
		ArrayList<String> final_hap_all_string_list=new ArrayList<String>();
		
		BufferedReader br_gcf = new BufferedReader(new FileReader(gcf_file));
		BufferedReader br_var = new BufferedReader(new FileReader(vars_intra_file));
		BufferedWriter bw_full_inter=new BufferedWriter(new FileWriter(all_inter_file));
		String curr_var_line = br_var.readLine();
		curr_var_line = br_var.readLine(); 
		while(curr_var_line!=null) {
			String[] curr_pos_arr = curr_var_line.split("\t");
			loci_string_list.add(curr_pos_arr[0]);
			ArrayList<String> new_var_list=new ArrayList<String>();
			complete_hap_seq_listlist.add(new_var_list); // num_var * num_hap
			curr_var_line = br_var.readLine();
		}
		br_var.close();
		String curr_gcf_line = br_gcf.readLine();
		while(curr_gcf_line!=null) {
			String[] hap_array = curr_gcf_line.split("\t");
			String[] every_position = hap_array[0].split("");
			ArrayList<String> tmp_hap_list=new ArrayList<String>();	
			String new_full_hap = "";
			for(int i=0;i<every_position.length;i++) {
				tmp_hap_list.add(every_position[i]);
				if(every_position[i].equals("0")||every_position[i].equals("1")) {
					new_full_hap=new_full_hap+every_position[i];
				}
			}
			hap_all_string_list.add(new_full_hap); 
			if(!tmp_hap_list.contains("?")) {
				complete_hap++;
				int var_pos=0;
				for(int i=0;i<tmp_hap_list.size();i++) {
					if(every_position[i].equals("0")||every_position[i].equals("1")) {
						complete_hap_seq_listlist.get(var_pos).add(every_position[i]);
						var_pos++;
					}
				}
			}
			curr_gcf_line = br_gcf.readLine();
		}
		//Write complete_inter_file
		if(complete_hap!=0) {
			BufferedWriter bw_complete_file=new BufferedWriter(new FileWriter(complete_inter_file));
			bw_complete_file.write("Hap_ID");
			for(int h=0;h<complete_hap_seq_listlist.get(0).size();h++) {
				bw_complete_file.write("\t"+"h"+h);
			}
			bw_complete_file.write("\n");
			bw_complete_file.write("Freq");
			for(int h=0;h<complete_hap_seq_listlist.get(0).size();h++) {
				bw_complete_file.write("\t"+0.0);
			}
			bw_complete_file.write("\n");
			for(int l=0;l<loci_string_list.size();l++) {
				bw_complete_file.write(loci_string_list.get(l));
				for(int h=0;h<complete_hap_seq_listlist.get(l).size();h++) {
					bw_complete_file.write("\t"+complete_hap_seq_listlist.get(l).get(h));
				}
				bw_complete_file.write("\n");
			}
			bw_complete_file.close();	
		}else {
			System.out.println("For "+gcf_file+", no complete haplotypes, no complete_inter_file is generate");
		}

		br_gcf.close();
		// get rid of the identical haplotype in the hap_full_string_list
		for(int h=0;h<hap_all_string_list.size();h++) {
			if(!final_hap_all_string_list.contains(hap_all_string_list.get(h))) {
				final_hap_all_string_list.add(hap_all_string_list.get(h));
			}
		}
		// Convert final_hap_full_string_list to full_hap_seq_listlist 
		for (int h=0; h < final_hap_all_string_list.size(); h++) {
			ArrayList<String> tmp_hap_string_list=new ArrayList<String>();
			for(int i=0; i < final_hap_all_string_list.get(h).split("").length; i++) {
				tmp_hap_string_list.add(final_hap_all_string_list.get(h).split("")[i]);
			}
			all_hap_seq_listlist.add(tmp_hap_string_list);
		}
		//Write full_inter_file
		bw_full_inter.write("Hap_ID");
		for(int h=0;h<final_hap_all_string_list.size();h++) {
			bw_full_inter.write("\t"+"h"+h);
		}
		bw_full_inter.write("\n");
		bw_full_inter.write("Freq");
		for(int h=0;h<final_hap_all_string_list.size();h++) {
			bw_full_inter.write("\t"+0.0);
		}
		bw_full_inter.write("\n");
		for(int l=0;l<loci_string_list.size();l++) {
			bw_full_inter.write(loci_string_list.get(l));
			for(int h=0;h<all_hap_seq_listlist.size();h++) {
				bw_full_inter.write("\t"+all_hap_seq_listlist.get(h).get(l));
			}
			bw_full_inter.write("\n");
		}
		bw_full_inter.close();
	}
	
	/* 
	 * For the non_perfect_data, one or more variant positions may not be called.
	 * But in order to compare the ori_inter_file with recon_inter_file, the
	 * number and position of variants must be the same.
	 * compare_loci() will compare the variants in the ori_inter_file and 
	 * recon_inter_file, if the recon_inter_file missed any of the variant positions
	 * or called any false_positive variant positions, a new version recon_inter_file 
	 * that contains all variant positions will overwrite the old one.        
	 */
	public static void compare_loci(String ori_inter_file, 
    		String recon_inter_file) throws IOException, InterruptedException{
    	 
    	 BufferedReader br_ori_inter = new BufferedReader(new FileReader(
    			 ori_inter_file));
    	 BufferedReader br_recon_inter = new BufferedReader(new FileReader(
    			 recon_inter_file));
    	 
    	 ArrayList<String> ori_variant_position_list = new ArrayList<>();
    	 ArrayList<String> recon_variant_position_list = new ArrayList<>();
    	 ArrayList<ArrayList<String>> hap_seq_listlist=new ArrayList<ArrayList<String>>();
		 ArrayList<String> hap_string_list=new ArrayList<String>();
		 
    	 //Read original_inter_file, and generate ori_variant_position_list
    	 String curr_ori_inter = br_ori_inter.readLine(); // read hap_ID line
    	 curr_ori_inter = br_ori_inter.readLine(); // read freq line
    	 curr_ori_inter = br_ori_inter.readLine(); // read third line
    	 while(curr_ori_inter != null) {
    		 String[] var_position_line = curr_ori_inter.split("\t");
    	 	 String[] var_position = var_position_line[0].split(";");
    	 	 ori_variant_position_list.add(var_position[1]);
    	 	 curr_ori_inter = br_ori_inter.readLine();
    	 }
    	 br_ori_inter.close();
    	 
    	 //Read reconstruct_inter_file
    	 String curr_recon_inter = br_recon_inter.readLine(); //read hap_ID line
    	 String[] hap_id_array = curr_recon_inter.split("\t");
		 for(int id=1; id < hap_id_array.length; id++) { 
			 ArrayList<String> new_hap_list=new ArrayList<String>();
			 hap_seq_listlist.add(new_hap_list);
		 }
    	 curr_recon_inter = br_recon_inter.readLine(); //read freq line
    	 String[] freq_array = curr_recon_inter.split("\t");
    	 double[] hap_freq_array = new double[freq_array.length-1];
		 for (int h =1; h < freq_array.length; h++ ) {
			 //hap_freq_list.add(Double.parseDouble(freq_array[h]));
			 hap_freq_array[h-1]=Double.parseDouble(freq_array[h]);
		 }
    	 curr_recon_inter = br_recon_inter.readLine(); //read the variant line
    	 int recon_num_loci = 0;
    	 while(curr_recon_inter != null) {
			 String[] var_position_line = curr_recon_inter.split("\t");
    	 	 String[] var_position = var_position_line[0].split(";");
    	 	 // get rid of false_positive variant positions
    	 	 if(ori_variant_position_list.contains(var_position[1])) {
    	 		   recon_num_loci++;
    	 		   recon_variant_position_list.add(var_position[1]);
    	 		   for (int index = 1; index < var_position_line.length; index ++) {
    	 			  hap_seq_listlist.get(index-1).add(var_position_line[index]);
    	 		   
    	 		   }
    	 	 }
    	 	 curr_recon_inter = br_recon_inter.readLine();
    	 }
    	 br_recon_inter.close();
    	 // change hap_seq_listlist to hap_string_list
		 for (int h=0; h<hap_seq_listlist.size();h++) {
			 String hap_str = "" ;
			 for (int index =0; index < hap_seq_listlist.get(h).size();index ++) {
				 hap_str = hap_str + hap_seq_listlist.get(h).get(index);
			 }
			 hap_string_list.add(hap_str);
		 }
		/*
		 * Check for identical haplotypes (Hap_ID), combine the haplotypes and
		 * corresponding frequency
		 */
		ArrayList<Integer> identical_hap= new ArrayList<>();
		ArrayList<String> final_hapString_list = new ArrayList<>();
		ArrayList<Double> final_freq_list = new ArrayList<>();
		for(int h1=0;h1<hap_string_list.size()-1;h1++) {
			for(int h2=h1+1;h2<hap_string_list.size();h2++) {
				if(hap_string_list.get(h2).equals(hap_string_list.get(h1))) {
					identical_hap.add(h2);
					hap_freq_array[h1]=hap_freq_array[h1]+hap_freq_array[h2];
				}
			}
		}
		for(int h=0;h<hap_string_list.size();h++) {
			 if(!identical_hap.contains(h)) {
				 final_hapString_list.add(hap_string_list.get(h));
				 final_freq_list.add(hap_freq_array[h]);
				 
			 }
		}
		/*
	     * Re_name the hap_ID according to haplotypes in the hap_string_list;
	     * The HapID is naturally coded by 0 and 1 strings. 
		 * For very long haplotypes, we hope to convert them in to an integer 
		 * of base 16. To indicate it is a haplotype, we add an "h" before the integer 
		 */
		for(int h=0;h<final_hapString_list.size();h++) {
			for(int locus=0;locus<recon_num_loci;locus++) {
				if(final_hapString_list.get(h).charAt(locus)!='1' && final_hapString_list.get(h).charAt(locus)!='0') {
					System.out.println("Exception: hap-ID "+final_hapString_list.get(h)+" is not 1/0 coded."
			+ "Returned without converting to base 16");
				}
			}
		}    
		// after check, convert the IDs to hexadecimal: 
		String[] new_hap_IDs = new String[final_hapString_list.size()];
		for(int h=0;h<final_hapString_list.size();h++) {
			BigInteger the_integer=new BigInteger(final_hapString_list.get(h),2);
		    new_hap_IDs[h]="h"+the_integer.toString(16);
		}
		
    	 //Over-write recon_inter_file, for those variant_positions 
    	 //in the ori_inter_file that are not called in the recon_inter_file, 
    	 //write 0 for all haplotypes; 
         PrintWriter pw_inter = new PrintWriter(new FileWriter(recon_inter_file, false)); 
         pw_inter.append("Hap_ID");
		 for (int h=0; h<new_hap_IDs.length;h++) {
			 pw_inter.append("\t"+new_hap_IDs[h]);
		 }
		 pw_inter.append("\n");
		 pw_inter.append("Freq");
		 for (int h=0; h<final_freq_list.size();h++) {
			 pw_inter.append("\t"+ final_freq_list.get(h));
		 }
		 pw_inter.append("\n");
		 int loci_position =0;
		 for (int l=0;l<ori_variant_position_list.size();l++) {
			 pw_inter.write("0;"+ori_variant_position_list.get(l)+";"
					 +ori_variant_position_list.get(l)+";0:1");
			 if(recon_variant_position_list.contains(ori_variant_position_list.get(l))) {
				 for(int h=0; h<final_hapString_list.size();h++) {
					 pw_inter.append("\t"+ final_hapString_list.get(h).charAt(loci_position));
				 }
				 loci_position++;
				 pw_inter.write("\n");
			 }else {
				 for(int h=0; h<final_hapString_list.size();h++) {
					 pw_inter.append("\t"+ "0");
				 }
				 pw_inter.write("\n");
			 }
		 }
		 pw_inter.close();
    }
	
	/*
	 * comparehapInter() and global_hap_evaluator():
	 * It's used to compare recon_inter file with ori_inter file.
	 * For the first graph_coloring results, it will only compare the haplotypes,
	 * won't compare the frequency. 
	 */
	public static void comparehapInter(String orig_inter_file, 
			String recon_inter_file, double quasi_cutoff, String output_files_prefix,
			String project_name) throws IOException, InterruptedException {
		double[] multi_pool_record = new double[5];
		multi_pool_record = global_hap_evaluator(orig_inter_file, 
				recon_inter_file, quasi_cutoff, output_files_prefix,
				project_name);
		PrintWriter pw2 = new PrintWriter(
	            new FileWriter(output_files_prefix + "_" + quasi_cutoff + 
	            		"_global_haps_average_results.txt", false));
	        pw2.append("## parameters: cut-off = " + quasi_cutoff + "\n");
	        pw2.append("#Project_Name\t"
	                + "Mean_OH_Recovery\t"
	                + "Mean_Mean_Min_Seq_Diff\t" // mean of means
	                + "Total_RH/Total_OH\t"
	                + "Total_OH\t"
	                + "Total_RH\n");
	            pw2.append(project_name + "\t");
	            for (int i = 0; i < 5; i++) {
	                pw2.append(multi_pool_record[i] + "\t");
	            }
	            pw2.append("\n");
	        pw2.close();
	}

	public static double[] global_hap_evaluator(String orig_inter_file, 
			String recon_inter_file, double quasi_cutoff, String output_files_prefix,
			String project_name) throws IOException, InterruptedException {
		
		HapConfigInterFile orig_haps = new HapConfigInterFile(orig_inter_file);
		HapConfigInterFile recon_haps = new HapConfigInterFile(recon_inter_file);
		int max_pos_diff = (int)Math.floor(orig_haps.num_loci*quasi_cutoff);
		double diff_ct = 0.0;
		int num_accurate = 0;
		HashSet<Integer> quasi_indices = new HashSet<Integer>();
		
		if (orig_haps.num_loci != recon_haps.num_loci) {
			System.out.println("Error: orig_haps.num_loci!=recon_haps.num_loci: \n"
	                + "Returned without comparison!");
		}
		ArrayList<String> ori_hap_id = new ArrayList<String>();
	       for (int ho = 0; ho < orig_haps.num_global_hap; ho++) {
	    	   ori_hap_id.add(orig_haps.hap_IDs[ho]);
	       }
		
		int[] min_diff_hap = new int[orig_haps.num_global_hap];
	    // Index of the closest reconstructed haplotype.
	    String[] min_diff_ID = new String[orig_haps.num_global_hap];
	    // Storing the IDs of matched haps in reconstructed HapConfig
	    int[] min_diff_pos = new int[orig_haps.num_global_hap];
	    // Number of differing loci to closest reconstruction.

        int[] max_diff_haps = new int[orig_haps.num_global_hap];
        int[][] pairwise_diff_num = new int[orig_haps.num_global_hap][recon_haps.num_global_hap];		
		
        for (int ori_h_index = 0; ori_h_index <orig_haps.num_global_hap; ori_h_index++) {
            int h_ori = orig_haps.hapID2index.get(orig_haps.hap_IDs[ori_h_index]);
            int min_hap_index = 0;
            int min_diff = orig_haps.num_loci;
            for (int h_recon = 0; h_recon < recon_haps.num_global_hap; h_recon++) {
                int diff = 0;
                for (int p = 0; p < orig_haps.num_loci; p++) {
                    if (Double.compare(orig_haps.global_haps[h_ori][p],
                        recon_haps.global_haps[h_recon][p]) != 0)
                        diff++;
                }
                pairwise_diff_num[ori_h_index][h_recon] = diff;
                if (diff < min_diff) {
                    min_hap_index = h_recon;
                    min_diff = diff;
                }
            } // end_of_for(h_recon)
            if (min_diff <= max_pos_diff)
                num_accurate++; 
            min_diff_hap[ori_h_index] = min_hap_index;
            min_diff_ID[ori_h_index] = recon_haps.hap_IDs[min_hap_index];
            min_diff_pos[ori_h_index] = min_diff;
            diff_ct += min_diff; // Cumulative min_diff
        }// end_of_for(ori_h_index)
        for (int h_id = 0; h_id < orig_haps.num_global_hap; h_id++) {
            for (int hr = 0; hr < recon_haps.num_global_hap; hr++) {
                if (pairwise_diff_num[h_id][hr] <= max_pos_diff) {
                    max_diff_haps[h_id]++;
                    quasi_indices.add(hr); 
                }
            }
        }
        PrintWriter pw1 = new PrintWriter(
            new FileWriter(output_files_prefix + "_" + quasi_cutoff + "_global_haplotypes.result.txt", false));
        pw1.append("Project_Name"+"\t"+"Ori_Hap_ID"+"\t"+"Closest_Recon_Hap_ID"+"\t"
            +"Min_diff_Pos"+"\t"+"Num_of_Recon_Meet_Cutoff\n");
        for (int h_index = 0; h_index < orig_haps.num_global_hap; h_index++) {
            pw1.append(project_name + "\t" + orig_haps.hap_IDs[h_index] + "\t"
                + min_diff_ID[h_index] + "\t" + min_diff_pos[h_index]
                + "\t" + max_diff_haps[h_index] + "\n");
        }
        pw1.close();

        return new double[] {
            (double) num_accurate / orig_haps.num_global_hap,
            diff_ct /orig_haps.num_global_hap,
            (double)recon_haps.num_global_hap/orig_haps.num_global_hap,
            orig_haps.num_global_hap,
            recon_haps.num_global_hap
        };
        
	}
	

	public static void main(String[] args)
			throws IOException, InterruptedException {
		
		String project_name= args[0];
    	double quasi_cutoff= Double.parseDouble(args[1]); 
    	int num_pools = Integer.parseInt(args[2]);
    	String gs_dir= args[3]+"/"; 
    	String output_dir= args[4]+"/"; 
    	String inter_dir = args[5] + "/";
    	String gcf_dir= inter_dir + "gcf/";
    	String vars_intra_file = inter_dir+project_name+"_vars.intra_freq.txt";
    	String all_inter_file;
    	String complete_inter_file;
    	String gcf_file;
    	String ori_inter_file= gs_dir+project_name+"_haps.inter_freq_vars.txt";	
    	
    	for(int p=0;p<num_pools;p++) {
    		//generate_inter_file()
    		gcf_file=gcf_dir+project_name+"_p"+p+".gcf";
        	all_inter_file = output_dir+project_name+"_p"+p+"_gc1_all.inter_freq_haps.txt";
        	complete_inter_file = output_dir+project_name+"_p"+p+"_gc1_complete.inter_freq_haps.txt";
        	generate_inter_file(gcf_file,vars_intra_file,all_inter_file,complete_inter_file);
        	// compare_loci()
        	String recon_inter_file = output_dir+project_name+"_p"+p+"_gc1_all.inter_freq_haps.txt";
        	compare_loci(ori_inter_file,recon_inter_file);
        	recon_inter_file = output_dir+project_name+"_p"+p+"_gc1_all.inter_freq_haps.txt";
        	String output_files_prefix = output_dir+project_name+"_p"+p+"_gc1_all";
        	comparehapInter(ori_inter_file,recon_inter_file, quasi_cutoff, 
        			output_files_prefix,project_name); 
        	File complete_file = new File(output_dir+project_name+"_p"+p+"_gc1_complete.inter_freq_haps.txt");
        	if(complete_file.exists()) {
        		recon_inter_file = output_dir+project_name+"_p"+p+"_gc1_complete.inter_freq_haps.txt";
            	compare_loci(ori_inter_file,recon_inter_file);
            	recon_inter_file = output_dir+project_name+"_p"+p+"_gc1_complete.inter_freq_haps.txt";
            	output_files_prefix = output_dir+project_name+"_p"+p+"_gc1_complete";
            	comparehapInter(ori_inter_file,recon_inter_file, quasi_cutoff, 
            			output_files_prefix,project_name); 
        	}
    	}
	}
}
