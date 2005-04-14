// Biclustering.java, copyright by Yizong Cheng, 2002
// Citation: Cheng Y and Church GM, "Biclustering of expression data",
// Proceedings of ISMB 2000, 93-103.
// Compile and run with Java 2
// Usage:  java Biclustering
// Record file name: Biclusters.txt
// Send questions and suggestions to cheng@uc.edu
// Available as http://cheng.ececs.uc.edu/biclustering/Biclustering.java
// Input file format: tab-delimited ASCII file as converted from Excel

//04-12-05 by yh: In two cases, the data is regarded as NA.
//1. blank on that column and row
//2. it's not a value, some characters, like 'NA'
//all comments embedded in the code are written by yh

import javax.swing.*;
import java.lang.*;
import java.util.*;
import java.io.*;
import java.awt.event.*;
import java.awt.*;


public class Biclustering extends JFrame implements MouseListener{

	int imagewidth = 1100;
	int imageheight = 800;
	JMenuBar menuBar = new JMenuBar();
	JMenu fileMenu = new JMenu("File");
	JMenuItem fileOpen = new JMenuItem("Open");
	JMenuItem clearRecord = new JMenuItem("Clear Record");
	JMenuItem fileExit = new JMenuItem("Exit");
	JMenu parameters = new JMenu("Parameters");
	JMenuItem score = new JMenuItem("Max H Score");
	JMenuItem minRows = new JMenuItem("Min Number of Rows");
	JMenuItem minCol = new JMenuItem("Min Number of Columns");
	JButton next = new JButton("Next Bicluster");
	JButton record = new JButton("Record");
	JPanel contents = null;
	Image offscreen = null;
	Graphics offGraphics = null;
	RepaintManager rm = null;
	String[] columnName = null;
	int numberOfColumns = 0;
	String[] rowName = null;
	int numberOfRows = 0;
	short[][] matrix = null;
	Random rand = new Random();
	int maxScore = 50;
	int minHeight = 5;
	int minWidth = 6;
	int batchThreshold = 100;
	boolean[] remainingR = null;
	boolean[] remainingC = null;
	double[] rowMean = null;
	double[] columnMean = null;
	double[] rowScore = null;
	double[] columnScore = null;
	double mean = 0;
	int smWidth = 0;
	int smHeight = 0;
	double HScore = 0;

	public Biclustering(){
		enableEvents(AWTEvent.WINDOW_EVENT_MASK);
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		setSize(screenSize);
		imagewidth = screenSize.width;
		imageheight = screenSize.height;
		setLocation(0,0);
		setTitle("Biclustering");

		fileExit.addActionListener(new ActionListener(){
			                           public void actionPerformed(ActionEvent e){
				                           System.exit(0);
			                           }
		                           });

		fileOpen.addActionListener(new ActionListener(){
			                           public void actionPerformed(ActionEvent e){
				                           openFile();
			                           }
		                           });

		clearRecord.addActionListener(new ActionListener(){
			                              public void actionPerformed(ActionEvent e){
				                              File file = new File("Biclusters.txt");
				                              if (file.exists()) file.delete();
			                              }
		                              });

		next.addActionListener(new ActionListener(){
			                       public void actionPerformed(ActionEvent e){
				                       getBicluster();
			                       }
		                       });

		record.addActionListener(new ActionListener(){
			                         public void actionPerformed(ActionEvent e){
				                         appendRecord();
			                         }
		                         });

		score.addActionListener(new ActionListener(){
			                        public void actionPerformed(ActionEvent e){
				                        boolean more = true;
				                        while (more){
					                        more = false;
					                        String inWord = JOptionPane.showInputDialog("Enter Max H Score");
					                        try{
						                        maxScore = Integer.parseInt(inWord);
					                        }catch(NumberFormatException e2){ more = true; }
				                        }
			                        }
		                        });

		minRows.addActionListener(new ActionListener(){
			                          public void actionPerformed(ActionEvent e){
				                          boolean more = true;
				                          while (more){
					                          more = false;
					                          String inWord = JOptionPane.showInputDialog("Enter Min Number of Rows");
					                          try{
						                          minHeight = Integer.parseInt(inWord);
					                          }catch(NumberFormatException e2){ more = true; }
				                          }
			                          }
		                          });

		minCol.addActionListener(new ActionListener(){
			                         public void actionPerformed(ActionEvent e){
				                         boolean more = true;
				                         while (more){
					                         more = false;
					                         String inWord = JOptionPane.showInputDialog("Enter Min Number of Columns");
					                         try{
						                         minWidth = Integer.parseInt(inWord);
					                         }catch(NumberFormatException e2){ more = true; }
				                         }
			                         }
		                         });

		fileMenu.add(fileOpen);
		fileMenu.add(clearRecord);
		fileMenu.add(fileExit);
		parameters.add(score);
		parameters.add(minRows);
		parameters.add(minCol);
		menuBar.add(fileMenu);
		menuBar.add(parameters);
		menuBar.add(next);
		menuBar.add(record);
		next.setVisible(false);
		record.setVisible(false);
		this.setJMenuBar(menuBar);
		contents = (JPanel)getContentPane();
		addMouseListener(this);
	}

	protected void processWindowEvent(WindowEvent e){
		super.processWindowEvent(e);
		if (e.getID() == WindowEvent.WINDOW_CLOSING)
			System.exit(0);
	}

	void openFile(){
		JFileChooser chooser = new JFileChooser();
		int result = chooser.showOpenDialog(this);
		int n = 0;
		int skipRows = 0;
		int skipColumns = 0;
		if (result == JFileChooser.APPROVE_OPTION){
			File file = chooser.getSelectedFile();
			BufferedReader br = null;
			try{
				br = new BufferedReader(new FileReader(file));
			}catch(FileNotFoundException e){
				JOptionPane.showMessageDialog(this, e.getMessage());
				return;
			}
			try{
				br.mark(1000);	//1000 is the readAheadLimit. approx the number of characters in one line.
				String line = br.readLine();  //one line has been read in already.
				StringTokenizer st = new StringTokenizer(line, "\t");
				n = st.countTokens();
				String inWord = null;
				boolean more = true;
				while (more){
					more = false;
					inWord = JOptionPane.showInputDialog("How many columns to skip?");
					try{
						skipColumns = Integer.parseInt(inWord);
					}catch(NumberFormatException e){ more = true; }
				}
				st = new StringTokenizer(line, "\t");
				for (int i = 0; i < skipColumns; i++) st.nextToken();
				numberOfColumns = n - skipColumns;
				columnName = new String[numberOfColumns];
				for (int i = 0; i < numberOfColumns; i++)
					columnName[i] = st.nextToken();
				more = true;
				while (more){
					more = false;
					inWord = JOptionPane.showInputDialog("How many rows to skip?");
					try{
						skipRows = Integer.parseInt(inWord);
					}catch(NumberFormatException e){ more = true; }
				}
				br.reset();	//bug here. reset the stream to the start of file. one line has been read before.
				for (int i = 1; i < skipRows; i++) br.readLine();
				numberOfRows = 0;
				while((line = br.readLine()) != null){
					numberOfRows++;
				}
				br.close();
			}catch(IOException e){
				JOptionPane.showMessageDialog(this, e.getMessage());
				return;
			}
			JOptionPane.showMessageDialog(this, "There are "
			                              + Integer.toString(numberOfRows) + " rows and "
			                              + Integer.toString(numberOfColumns) + " columns.");
			try{
				br = new BufferedReader(new FileReader(file));
			}catch(FileNotFoundException e){
				JOptionPane.showMessageDialog(this, e.getMessage());
				return;
			}
			matrix = new short[numberOfRows][numberOfColumns];

			try{
				for (int i = 0; i < skipRows; i++) br.readLine();
				for (int i = 0; i < numberOfRows; i++){
					String line = br.readLine();
					String w = null;
					int tab = 0;
					for (int j = 0; j < skipColumns; j++)
						tab = line.indexOf('\t', tab) + 1;
					for (int j = 0; j < numberOfColumns; j++){
						int nextTab = line.indexOf('\t', tab);
						if (nextTab == tab + 1)
							matrix[i][j] = (short)(rand.nextInt(1600) - 800);	//04-12-05 by yh: regarded as NA(random)
						else{
							if (nextTab == -1) w = line.substring(tab);
							else w = line.substring(tab, nextTab);
							try{
								double d = Double.parseDouble(w);
								matrix[i][j] = (short)(d * 100.0);
							}catch(NumberFormatException e){
								matrix[i][j] = (short)(rand.nextInt(1600) - 800);	//04-12-05 by yh: regarded as NA(random)
							}
						}
						tab = nextTab + 1;
					}

				}
				br.close();
			}catch(IOException e){
				JOptionPane.showMessageDialog(this, e.getMessage());
				return;
			}
		}
		next.setVisible(true);
		remainingR = new boolean[numberOfRows];
		remainingC = new boolean[numberOfColumns];
		rowMean = new double[numberOfRows];
		columnMean = new double[numberOfColumns];
		rowScore = new double[numberOfRows];
		columnScore = new double[numberOfColumns];
		contents.setDoubleBuffered(true);
		rm = RepaintManager.currentManager(contents);
		offscreen = rm.getOffscreenBuffer(contents, imagewidth, imageheight);
		offGraphics = offscreen.getGraphics();
		
		//04-13-05 output the matrix
		outputMatrix();

	}

	public void scoring(){
		mean = 0;
		for (int j = 0; j < numberOfColumns; j++) if (remainingC[j])
				columnMean[j] = 0;	//initialization
		for (int i = 0; i < numberOfRows; i++) if (remainingR[i]){
				rowMean[i] = 0;	//initialization
				for (int j = 0; j < numberOfColumns; j++) if (remainingC[j]){
						rowMean[i] += matrix[i][j];
						columnMean[j] += matrix[i][j];
					}
				mean += rowMean[i];
				rowMean[i] /= smWidth;	//smWidth determined in getBicluster()
			}
		for (int j = 0; j < numberOfColumns; j++) if (remainingC[j])
				columnMean[j] /= smHeight;	//smHeight determined in getBicluster()
		mean /= smWidth * smHeight;
		HScore = 0;
		for (int j = 0; j < numberOfColumns; j++) if (remainingC[j])
				columnScore[j] = 0;	//initialization
		for (int i = 0; i < numberOfRows; i++) if (remainingR[i]){
				rowScore[i] = 0;	//initialization
				for (int j = 0; j < numberOfColumns; j++) if (remainingC[j]) {
						double r = matrix[i][j] - rowMean[i] - columnMean[j] + mean;
						r = r * r;
						rowScore[i] += r;
						columnScore[j] += r;
					}
				HScore += rowScore[i];
				rowScore[i] /= smWidth;
			}
		HScore /= smWidth * smHeight;
		for (int j = 0; j < numberOfColumns; j++) if (remainingC[j])
				columnScore[j] /= smHeight;
	}

	public void getBicluster(){
		for (int i = 0; i < numberOfRows; i++) remainingR[i] = true;
		for (int j = 0; j < numberOfColumns; j++) remainingC[j] = true;
		smWidth = numberOfColumns;
		smHeight = numberOfRows;
		scoring();
		int index = 0;
		while ((HScore > maxScore) && (index > -1)){	//if no more columns or rows can be removed,
			if (smHeight > batchThreshold){			//This algorithm 2 of the paper. Multiple node deletion, but kind of different,
											//no batchThreshold in the paper's algorithm
				for (int i = 0; i < numberOfRows; i++)
					if (remainingR[i] && (rowScore[i] > HScore)){
						remainingR[i] = false;
						smHeight--;
					}
			}else{	//find the maximum from rowScore and columnScore
				double ms = 0;
				index = -1;
				boolean row = true;
				if (smHeight > minHeight){
					for (int i = 0; i < numberOfRows; i++)
						if (remainingR[i] && (rowScore[i] > ms)){
							ms = rowScore[i];
							index = i;
						}
				}
				if (smWidth > minWidth){
					for (int i = 0; i < numberOfColumns; i++)
						if (remainingC[i] && (columnScore[i] > ms)){
							ms = columnScore[i];
							index = i;
							row = false;
						}
				}
				if (index > -1)
					if (row){
						remainingR[index] = false;
						smHeight--;
					}else{
						remainingC[index] = false;
						smWidth--;
					}
			}
			scoring();
		}

		display();
		//04-12-05 by yh: put random numbers in the block(cluster)
		for (int i = 0; i < numberOfRows; i++) if (remainingR[i])
				for (int j = 0; j < numberOfColumns; j++) if (remainingC[j])
						matrix[i][j] = (short)(rand.nextInt(1600) - 800);
		if (!record.isVisible()) record.setVisible(true);
	}

	public void display(){
		int[] sorted = new int[smWidth];
		int index = 0;
		for (int j = 0; j < numberOfColumns; j++) if (remainingC[j]) sorted[index++] = j;

		offGraphics.setColor(new Color(204, 204, 204));
		offGraphics.fillRect(0,0,imagewidth, imageheight);
		offGraphics.setFont(new Font("Arial", Font.PLAIN, 10));
		offGraphics.setColor(Color.black);
		int stepSize = imagewidth / (smWidth + 3);
		int nameOrder = 0;
		for (int i = 0; i < numberOfRows; i++) if (remainingR[i]){
				offGraphics.drawRect(stepSize, 400 - matrix[i][sorted[0]], 3, 3);
				for (int j = 1; j < smWidth; j++){
					offGraphics.drawRect(stepSize * (j + 1), 400 - matrix[i][sorted[j]], 3, 3);
					offGraphics.drawLine(stepSize * j, 400 - matrix[i][sorted[j - 1]],
					                     stepSize * (j + 1), 400 - matrix[i][sorted[j]]);
				}
			}
		offGraphics.drawString("H Score = " + Double.toString(HScore), 100, 30);	//04-14-05, output the Hscore in double, not integer
		Graphics g = contents.getGraphics();
		g.clipRect(0,0,imagewidth, imageheight);
		g.drawImage(offscreen, 0, 0, contents);
	}

	public void appendRecord(){
		PrintWriter recordFile = null;
		try{
			recordFile = new PrintWriter(new FileWriter("Biclusters.txt", true));
		}catch(IOException e){
			System.err.println(e.getMessage());
			record.setVisible(false);
			return;
		}
		
		recordFile.println(smHeight + " " + smWidth + " " + HScore);	//output the score in double
		for (int i = 0; i < numberOfRows; i++) if (remainingR[i])
				recordFile.print(i + " ");
		recordFile.println();
		for (int j = 0; j < numberOfColumns; j++) if (remainingC[j])
				recordFile.print(j + " ");
		
		recordFile.println();
		recordFile.close();
	}

	public void outputMatrix(){
	/*
	*04-13-05
	*	output the matrix that's read in, for comparison with biclustering.cc
	*/
		PrintWriter recordFile = null;
		try{
			recordFile = new PrintWriter(new FileWriter("Biclusters.txt", true));
		}catch(IOException e){
			System.err.println(e.getMessage());
			record.setVisible(false);
			return;
		}
		
		for (int i = 0; i < numberOfRows; i++)
		{
			for (int j = 0; j < numberOfColumns; j++)
				recordFile.print(matrix[i][j]+"\t");
			recordFile.println();
			}
		recordFile.println();
		recordFile.close();
	}

public void mouseClicked(MouseEvent e){}
	public void mouseEntered(MouseEvent e){}
	public void mouseExited(MouseEvent e){}
	public void mousePressed(MouseEvent e){}
	public void mouseReleased(MouseEvent e){}

	public static void main(String[] args){
		Biclustering bc = new Biclustering();
		bc.setVisible(true);
	}
}
