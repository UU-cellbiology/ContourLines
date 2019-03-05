package fiji.plugin.ContourLines;



public class Point implements Cloneable {
	/**X and Y coordinate of a point **/
	public float[] coords;
	
	//default constructor
	public Point()
	{
		coords = new float[2];
	}
	public Point(float x, float y)
	{
		coords = new float[2];
		coords[0]=x;
		coords[1]=y;
	}

	public Point clone()
	{
		
		Point pointret = new Point(coords[0],coords[1]);
		return pointret;
	}
	
	public double getX() {
		// TODO Auto-generated method stub
		return coords[0];
	}
	
	public double getY() {
		// TODO Auto-generated method stub
		return coords[1];
	}
	public float distance(float x, float y)
	{
		return (float) Math.sqrt(Math.pow(x-coords[0], 2)+Math.pow(y-coords[1], 2));
	}
}
