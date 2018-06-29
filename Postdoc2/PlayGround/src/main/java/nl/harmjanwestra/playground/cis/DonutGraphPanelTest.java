package nl.harmjanwestra.playground.cis;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.themes.DefaultTheme;
import nl.harmjanwestra.utilities.graphics.themes.Theme;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.awt.*;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Locale;
import java.util.regex.Pattern;

public class DonutGraphPanelTest {
	
	
	public static void main(String[] args) {
		
		Grid g = new Grid(200, 200, 7, 7, 20, 20);
		ArrayList<String> abbs = new ArrayList<String>();
		TextFile tf = null;
		
		try {
			tf = new TextFile("D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-06-05-cis-gtexrepl\\wofdr\\concordanceOfTopFx-nonfdr.txt", TextFile.R);
			tf.readLine();
			
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				DonutGraphPanel p = new DonutGraphPanel(1, 1);
				
				String abb = elems[1];
				Double all = Double.parseDouble(elems[3]);
				double overlap = Double.parseDouble(elems[8]);
				double concordant = Double.parseDouble(elems[15]);
				
				double perc0 = overlap / all;
				double perc1 = concordant / all;
				
				p.setData(new double[]{perc0, perc1});
				
				String n = NumberFormat.getNumberInstance(Locale.US).format(all);
				p.setText(abb + "\n(n=" + n + ")");
				p.setTheme(new eQTLGenTheme());
				abbs.add(elems[1] + ": " + elems[0].replaceAll("_", " ").replaceAll(".v7.egenes.txt.gz", ""));
				g.addPanel(p);
				elems = tf.readLineElems(TextFile.tab);
			}
			
			
			tf.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		DefaultTheme t = new DefaultTheme();
		
		try {
			g.draw("D:\\TMP\\test.pdf");
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
		
		System.out.println(Strings.concat(abbs, Pattern.compile(", ")));
	}
	
	static class eQTLGenTheme implements Theme {
		
		public final Font LARGE_FONT = new Font("Helvetica", 0, 14);
		public final Font LARGE_FONT_BOLD = new Font("Helvetica", 1, 14);
		public final Font MEDIUM_FONT = new Font("Helvetica", 0, 12);
		public final Font MEDIUM_FONT_BOLD = new Font("Helvetica", 1, 12);
		public final Font SMALL_FONT = new Font("Helvetica", 0, 10);
		public final Font SMALL_FONT_BOLD = new Font("Helvetica", 1, 10);
		public final Stroke strokeDashed = new BasicStroke(1.0F, 1, 1, 0.0F, new float[]{4.0F}, 0.0F);
		public final Stroke stroke2pt = new BasicStroke(2.0F, 1, 1);
		public final Stroke stroke2ptDashed = new BasicStroke(2.0F, 1, 1, 0.0F, new float[]{4.0F}, 0.0F);
		public final Stroke stroke = new BasicStroke(1.0F, 1, 1);
		private final Color darkgrey = new Color(100, 100, 100);
		private final Color lightgrey = new Color(225, 225, 225, 75);
		private final Color[] colors = new Color[]{
				new Color(255, 0, 0),
				new Color(0, 0, 255),
				new Color(98, 182, 177),
				new Color(116, 156, 80),
				new Color(124, 87, 147),
				new Color(174, 164, 140),
				new Color(185, 113, 65),
				new Color(63, 93, 126),
				new Color(79, 96, 64),
				new Color(109, 54, 96),
				new Color(110, 36, 30), new Color(32, 79, 74), new Color(81, 94, 28), new Color(67, 42, 82)};
		
		public eQTLGenTheme() {
		}
		
		public Color getColor(int i) {
			return this.colors[i % this.colors.length];
		}
		
		public Color getLightGrey() {
			return this.lightgrey;
		}
		
		public Color getDarkGrey() {
			return this.darkgrey;
		}
		
		public Font getLargeFont() {
			return this.LARGE_FONT;
		}
		
		public Font getMediumFont() {
			return this.MEDIUM_FONT;
		}
		
		public Font getLargeFontBold() {
			return this.LARGE_FONT_BOLD;
		}
		
		public Font getMediumFontBold() {
			return this.MEDIUM_FONT_BOLD;
		}
		
		public Font getSmallFont() {
			return this.SMALL_FONT;
		}
		
		public Font getSmallFontBold() {
			return this.SMALL_FONT_BOLD;
		}
		
		public Stroke getStroke() {
			return this.stroke;
		}
		
		public Stroke getStrokeDashed() {
			return this.strokeDashed;
		}
		
		public Stroke getThickStroke() {
			return this.stroke2pt;
		}
		
		public Stroke getThickStrokeDashed() {
			return this.stroke2ptDashed;
		}
		
		public Color getColorSetOpacity(int i, float v) {
			Color c = this.colors[i];
			int r = c.getRed();
			int g = c.getGreen();
			int b = c.getBlue();
			int a = (int) Math.floor((double) (v * 255.0F));
			return new Color(r, g, b, a);
		}
		
		public Color getDarkerColor(Color color, double perc) {
			double delta = 1.0D - perc;
			int r = (int) Math.ceil((double) color.getRed() * delta);
			int g = (int) Math.ceil((double) color.getGreen() * delta);
			int b = (int) Math.ceil((double) color.getBlue() * delta);
			return new Color(r, g, b);
		}
	}
}
