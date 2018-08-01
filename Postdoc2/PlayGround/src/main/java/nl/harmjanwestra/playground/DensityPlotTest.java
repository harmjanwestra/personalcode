package nl.harmjanwestra.playground;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.playground.cis.DensityGraphPanel;
import nl.harmjanwestra.utilities.graphics.Grid;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.IOException;
import java.util.ArrayList;

public class DensityPlotTest {
	
	
	public static void main(String[] args) {
		
		ArrayList<Double> x = new ArrayList<Double>();
		ArrayList<Double> y = new ArrayList<Double>();
		
		NormalDistribution dist = new NormalDistribution(0, 1);
		
		for (int i = 0; i < 100000; i++) {
			
			double val = dist.sample();
			
			double valy = dist.sample();
			x.add(val);
			y.add(valy);
		}
		
		
		Grid g = new Grid(250, 250, 1, 1, 50, 50);
		DensityGraphPanel p = new DensityGraphPanel(1, 1);
		p.setData(x, y);
		
		p.setPlotElems(true, false);
		g.addPanel(p);
		
		
		try {
			g.draw("D:\\TMP\\density.png");
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
		
	}
	
}

