import java.io.IOException;

public class OpenFPM3D {
	
    public OpenFPM3D() {
		System.out.println( "Loading OpenFPM3D" );
		    try {
		    	// Check native system
		        System.loadLibrary("openfpm3D"); // used for tests. This library in classpath only
		    } catch (UnsatisfiedLinkError e) {
		        try {
		            NativeUtils.loadLibraryFromJar("/libopenfpm3D.jnilib"); // during runtime. .DLL within .JAR
		        } catch (IOException e1) {
		            throw new RuntimeException(e1);
		        }
		    }

		System.out.println( "OpenFPM3D" );
	}
    
    public static native String echo(String in);

    public static native void setDT( double dt );
    public static native void setNeighborhoodRadius( double r );
    public static native void setBoundaryWidth( double w );
    public static native void setBoundaryHeight( double h );
    public static native void setBoundaryDepth( double d );

    // True = periodic, false = nonperiodic
    public static native void setBoundaryConditionX( boolean boundaryCondition );
    public static native void setBoundaryConditionY( boolean boundaryCondition );
    public static native void setBoundaryConditionZ( boolean boundaryCondition );

    // In these simulations we only use OpenFPM to manage positions and position-related mechanisms
    // A particle then contains a Java Object that stores simulation specific state per particle (like velocity, force, etc.)
    
    public static native void init( );
    public static native void initBoundary( );    
    public static native int addParticle( double[] position, Object particleState );

    public static native void map( );    // This name really bothers me

    public static native void updateCellList( );

    public static native double[] getParticlePosition( int key );
    public static native void setParticlePosition( int key, double[] newPosition );    
    public static native Object getParticleState( int key );
    public static native int[] getParticleNeighbors( int key );    
    
    // Don't forget about checking if neighbor==self
    // DESTROY

}
