package MiscFunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import PoolHap.LocusAnnotation;

public class HapConfigInterFile {
	    public int num_global_hap; // number of global haplotypes in this region
	    public int num_loci; // number of loci in this region

	    // Map the pool IDs to their indices in this.hap_IDs.
	    public HashMap<String, Integer> hapID2index;
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

	
	public HapConfigInterFile(String global_hap_input_file) throws IOException,
														InterruptedException {
      // Parse global_hap_input_file.
      BufferedReader br = new BufferedReader(new FileReader(global_hap_input_file));
      String[] header_ids = br.readLine().split("\t");
      String[] header_freqs = br.readLine().split("\t");
      this.num_global_hap = header_ids.length - 1;
      this.hap_IDs = new String[this.num_global_hap];
      this.global_haps_freq = new double[this.num_global_hap];
      this.hapID2index = new HashMap<String, Integer>();
      for (int h = 0; h < num_global_hap; h++) {
          this.hap_IDs[h] = header_ids[h + 1];
          this.hapID2index.put(this.hap_IDs[h], h);
          this.global_haps_freq[h] = Double.parseDouble(header_freqs[h + 1]);
      }    

      String line = br.readLine();
      while (line != null) {
          this.num_loci++;
          line = br.readLine();
      }
      br.close();
      this.global_haps_string = new String[this.num_global_hap][this.num_loci];
      this.global_haps = new double[this.num_global_hap][this.num_loci];
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
              this.global_haps[h_index][loci_index]=Double.parseDouble(tmp[h_index+1]);
          }
          loci_index++;
          line = br.readLine();
      }
      br.close();
	}
}
