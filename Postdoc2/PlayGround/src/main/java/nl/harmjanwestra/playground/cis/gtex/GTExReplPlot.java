package nl.harmjanwestra.playground.cis.gtex;

import umcg.genetica.graphics.panels.Panel;
import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;

import java.awt.*;

public class GTExReplPlot {
	
	public static void main(String[] args) {
	
	}
	
	public void run() {
	
	}
	
	public class StackedBarPlotPanel extends Panel {
		
		private String[] stacklabels; // [cats][stacks]
		private double[][][] stacks; // [cats][stacks][values]
		private String[] categoryLabels;
		private Range datarange;
		private boolean additive = false; // consider stacks to be cumulative (true), or to be subsets of each other (false)
		private boolean scalepercategory = false; // scale the values to 100% per category, or rather let the max value over all categroies be 100%
		
		public StackedBarPlotPanel(int nrRows, int nrCols) {
			super(nrRows, nrCols);
		}
		
		public void setStacklabels(String[] labels, int nrcategories) {
			this.stacklabels = labels;
		}
		
		public void setCategoryLabels(String[] labels) {
			this.categoryLabels = labels;
		}
		
		
		public void setData(double[][][] data) {
			this.stacks = data;
		}
		
		@Override
		public void draw(DefaultGraphics defaultGraphics) {
			
			// max X
			int nrcats = stacks.length;
			
			double[] maxys = new double[stacks.length];
			
			// determine max Y
			int catctr = 0;
			for (double[][] cat : stacks) {
				double maxy = 0;
				double sum = 0;
				if (cat != null) {
					for (double[] stack : cat) {
						if (stack != null) {
							for (double val : stack) {
								sum+=val;
								if (val > maxy) {
									maxy = val;
								}
							}
						}
					}
				}
				
				if(additive){
					maxys[catctr] = sum;
				} else {
					maxys[catctr] = maxy;
				}
				
				
				catctr++;
			}
			
			if (!scalepercategory) {
				double max = 0;
				for (double d : maxys) {
					if (d > max) {
						max = d;
					}
				}
				for (int q = 0; q < maxys.length; q++) {
					maxys[q] = max;
				}
				
				// scaling over all categories sets the max at the highest value
				datarange = new Range(0, 0, 1, max);
			} else {
				// scaling per category sets the max at 100%
				datarange = new Range(0, 0, 1, 1);
			}
			
			
			
			Graphics2D g2d = defaultGraphics.getG2d();
			
			
			
			
			
		}
	}
	
}
