package MiscFunctions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;


public class InterRegionsSplit {
	
	/*
	 * generate_split_interfile():
	 * Read ori_inter_file, and only read those information within this region;
	 * Combine identical haplotypes and their corresponding frequency within 
	 * this region; Write a new regional_inter_file under gold_standard.
	 * 
	 * TODO: For now, this can only work for perfect data, when the variant positions
	 * in the "_vars.intra_freq" are the same as in the ori_inter_file
	 * 
	 */
	 public static void generate_split_interfile(
		        String project_name,
		        String gs_dir,
		        String curr_region,
		        int level,
		        int region_count) throws IOException, InterruptedException {
		 
		 ArrayList<ArrayList<String>> hap_seq_listlist=new ArrayList<ArrayList<String>>();
		 ArrayList<String> hap_string_list=new ArrayList<String>();
		 ArrayList<String> final_hap_string_list=new ArrayList<String>();
		 ArrayList<ArrayList<String>> final_hap_seq_listlist=new ArrayList<ArrayList<String>>();
		 ArrayList<String> loci_string_list=new ArrayList<String>();
		 ArrayList<Double> hap_freq_list=new ArrayList<Double>();
		 HashMap<String, Double> hap2fre = new HashMap<String, Double>();
		 
		 String gs_inter_file = gs_dir + project_name + "_haps.inter_freq_vars.txt";
		 String region_inter_file = gs_dir + project_name + "_level_" + level + 
				 "_region_" + region_count + "_haps.inter_freq_vars.txt";
		 String[] region_arr = curr_region.split(":");
		 int loci_begin = Integer.parseInt(region_arr[0]);
		 int loci_end = Integer.parseInt(region_arr[1]);
		 BufferedReader br_inter = new BufferedReader(new FileReader(gs_inter_file));
		 String curr_line = br_inter.readLine(); // read the header, hap_ID
		 String[] hap_id_array = curr_line.split("\t");
		 for(int id=1; id < hap_id_array.length; id++) { 
			 ArrayList<String> new_hap_list=new ArrayList<String>();
			 hap_seq_listlist.add(new_hap_list);
		 }
		 curr_line = br_inter.readLine(); // read the freq line
		 String[] freq_array = curr_line.split("\t");
		 for (int h =1; h < freq_array.length; h++ ) {
			 hap_freq_list.add(Double.parseDouble(freq_array[h]));
		 }
		 for (int l =0; l < loci_begin; l++) { //read those lines before the loci_begin
			 curr_line = br_inter.readLine();
		 }
		 for (int l=0; l <= (loci_end - loci_begin); l++) { // read those lines from the loci_begin to loci_end
			 curr_line = br_inter.readLine();
			 String[] loci_array = curr_line.split("\t");
			 loci_string_list.add(loci_array[0]);
			 for (int index = 1; index < loci_array.length; index ++) {
				 hap_seq_listlist.get(index-1).add(loci_array[index]);
			 } 
		 }
		 br_inter.close();
		 // change hap_seq_listlist to hap_string_list
		 for (int h=0; h<hap_seq_listlist.size();h++) {
			 String hap_str = "" ;
			 for (int index =0; index < hap_seq_listlist.get(h).size();index ++) {
				 hap_str = hap_str + hap_seq_listlist.get(h).get(index);
			 }
			 hap_string_list.add(hap_str);
		 }
		 // generate hap2fre hashmap
		 // Combine identical haplotypes and their corresponding frequency
		 for (int h=0; h<hap_string_list.size();h++) {
			 if(!hap2fre.containsKey(hap_string_list.get(h))) {
				hap2fre.put(hap_string_list.get(h), hap_freq_list.get(h));
			}else if (hap2fre.containsKey(hap_string_list.get(h))){
				hap2fre.put(hap_string_list.get(h), 
						(hap2fre.get(hap_string_list.get(h))+hap_freq_list.get(h)));
			}
		 }
		 // generate final_hap_string_list using hap2fre
		for ( String key : hap2fre.keySet() ) {
		    final_hap_string_list.add(key);
		}
		// final_hap_string_list transfer to final_hap_seq_listlist
		// final_hap_seq_listlist contains each haplotypes and their loci(0/1)
		for (int h=0; h < final_hap_string_list.size(); h++) {
			ArrayList<String> tmp_hap_string_list=new ArrayList<String>();
			for(int i=0; i < final_hap_string_list.get(h).split("").length; i++) {
				tmp_hap_string_list.add(final_hap_string_list.get(h).split("")[i]);
			}
			final_hap_seq_listlist.add(tmp_hap_string_list);
		}
		 // write output inter file for each regions
		 BufferedWriter bw_inter=new BufferedWriter(new FileWriter(
				 region_inter_file));
		 bw_inter.write("Hap_ID");
		 for (int h=0; h<final_hap_string_list.size();h++) {
			 bw_inter.write("\t"+"h"+h);
		 }
		 bw_inter.write("\n");
		 bw_inter.write("Freq");
		 for (int h=0; h<final_hap_string_list.size();h++) {
			 double curr_hap_freq = hap2fre.get(final_hap_string_list.get(h));
			 bw_inter.write("\t"+ curr_hap_freq);
		 }
		 bw_inter.write("\n");
		 for (int l=0;l<loci_string_list.size();l++) {
			 bw_inter.write(loci_string_list.get(l));
			 for(int h=0; h<final_hap_seq_listlist.size();h++) {
				 bw_inter.write("\t"+ final_hap_seq_listlist.get(h).get(l));
			 }
			 bw_inter.write("\n");
		 }
		 bw_inter.close();
	 }
	 
	public static void main(String[] args) throws IOException, InterruptedException{
		String project_name= args[0];
    	String gs_dir= args[1]+"/";
    	String inter_dir= args[2]+"/";
    	String dc_plan_file=inter_dir+project_name+"_dc_plan.txt";
        BufferedReader br_dc = new BufferedReader(new FileReader(dc_plan_file));
        String curr_line = br_dc.readLine();//read first line "Level I"
        curr_line = br_dc.readLine(); // read the regions for level I
        String[] regions_arr = curr_line.split("\t");
        for(int r=0; r<regions_arr.length;r++) {
        	String curr_region = regions_arr[r]; // for exa: 0:6
        	int level = 1;
        	int region_count = r;
        	generate_split_interfile(project_name,gs_dir,curr_region,level,region_count);
        }
        curr_line = br_dc.readLine();//read first line "Level II"
        curr_line = br_dc.readLine(); // read the regions for level II
        regions_arr = curr_line.split("\t");
        for(int r=0; r<regions_arr.length;r++) {
        	String curr_region = regions_arr[r]; // for exa: 3:9
        	int level = 2;
        	int region_count = r;
        	generate_split_interfile(project_name,gs_dir,curr_region,level,region_count);
        }
        br_dc.close();
	}

}
