public class Main {
	
	public static void main(String[] args) {
		String[] months = new String[]{
				"Dec",
				"Jan",
				"Feb",
				"Mrt",
				"Apr",
				"Mar",
				"Jun",
				"Jul",
				"Aug",
				"Sep",
				"Oct",
				"Nov",
				"Dec",
				"Jan",
				"Feb",
				"Mrt",
				"Apr",
				"Mar",
				"Jun",
				"Jul",
				"Aug",
				"Sep",
				"Oct",
				"Nov",
				"Dec"
		};
		
		int[] extrainkomsten = new int[months.length];
		extrainkomsten[0] = 220;
		extrainkomsten[4] = 3513 + 500; // saldo amerika + belastingteruggaaf (??)
		
		
		String[] names = new String[]{
				"Salaris",
				"Huur",
				"Trein",
				"Eten",
				"Internet",
				"Verzekering",
				"Telefoon",
				"Duo"
		};
		
		int[] kosten = new int[]{
				2561,
				-770,
				-317,
				-600,
				-35,
				-110,
				-30,
				-220
		};
		
		int delta = 0;
		for (int q = 0; q < kosten.length; q++) {
			delta += kosten[q];
		}
		
		int prevsaldo = 250;
		System.out.println("Standaard delta: " + delta);
		for (int maand = 0; maand < months.length; maand++) {
			int newsaldo = prevsaldo + delta + extrainkomsten[maand];
			System.out.println(months[maand] + "\t" + prevsaldo + "\t" + newsaldo);
			prevsaldo = newsaldo;
		}
	}
}
