package Viral_Reconstructions_Tools;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.NumberFormat;


public class haplofile_nloci_2nhaps {
	
	public static void main(String[] args) {
		//String folder="C:\\Users\\lenovo\\Desktop\\HIV project\\Hippo\\";
		String output_haplofile="haplo_file_n_"+args[0]; // "num_loci"
		int num_loci=Integer.parseInt(args[0]);
		haplo_file(num_loci,output_haplofile);
	}
	public static void haplo_file(int num_of_loci,String output_haplofile) {
		try {
			BufferedWriter bw_file=new BufferedWriter(new FileWriter(output_haplofile));	
			double haps_2n=1.0; // Because use int will be out of range, I use double instead
			for(int i=0;i<num_of_loci;i++) {
				haps_2n=haps_2n*2;
			}
			// To be notice: By using NumberFormat, the double won't show as 1.2099E12 
			NumberFormat nf = NumberFormat.getInstance();
			nf.setGroupingUsed(false);
			System.out.println(nf.format(haps_2n));
			bw_file.write(nf.format(haps_2n)+"\n");
			for(int hap=0;hap<haps_2n;hap++) {
				String vc_str = Integer.toBinaryString(hap);
				String[] vc_arr= vc_str.split("");
				for(int num_0=0;num_0<num_of_loci - vc_arr.length;num_0++) {
					bw_file.write("0"+" ");
				}
				for(int loci=0;loci<vc_arr.length-1;loci++) {
					bw_file.write(vc_arr[loci]+" ");
				}
				bw_file.write(vc_arr[vc_arr.length-1]+"\n");
			}
			bw_file.close();
		}catch(Exception e) {e.printStackTrace();}
	}
}