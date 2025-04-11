"use strict";

{
    const { chroma, qx, require, THREE, async, desk } = window;

    const ACVD = {};

    ACVD.extractMeshes = async function ( file, objects, viewer ) {

        const meshes = {};

        const promises = Object.keys( objects ).map( async function ( name ) {

            let opts = JSON.parse( JSON.stringify(objects[ name ]) );
            opts.action = "extract_meshes";
            opts.input_volume = file;
            opts.max_number_of_vertices = opts.numberOfVertices || opts.vertices || 10000;

            const res = await desk.Actions.executeAsync( opts );
            const meshFile = res.outputDirectory + '0.vtk';
            meshes[ name ] = { file : meshFile };
            if ( !viewer ) return;

            try {

                const mesh = await viewer.addFileAsync( meshFile, { label : name } );
                if ( opts.material ) mesh.material.setValues( opts.material );
                if ( opts.transparent != undefined ) mesh.transparent = opts.transparent;
                if ( opts.renderOrder != undefined ) mesh.renderOrder = opts.renderOrder;
                meshes[ name ].mesh = mesh;

            } catch ( e ) { }

        } );

        await Promise.all( promises );
        return meshes;

    } ;

    ACVD.cleanMeshes = async function ( inputDir, outputDir ) {

        const path = require( "path" );
        await desk.FileSystem.traverseAsync( inputDir, async function ( inputMesh, cb ) {

            const clean = await desk.Actions.executeAsync( {
                action : "mesh2ply",
                inputMesh,
                clean : true
            } );

            const relative = path.relative( inputDir, inputMesh );
            const destination = path.join( outputDir, relative );
            const destDir = path.dirname( destination );
            await desk.FileSystem.mkdirpAsync( destDir );

            await desk.Actions.executeAsync( {
                action : "copy",
                source : path.join( clean.outputDirectory, "mesh.ply" ),
                destination
            } );

            cb();

        }, true );

    };

    ACVD.batchMeshConvert = async function ( inputDir, outputDir, format = "stl" ) {

        const path = require( "path" );
        await desk.FileSystem.traverseAsync( inputDir, async function ( inputMesh, cb ) {

            const clean = await desk.Actions.executeAsync( {
                action : "mesh2" + format,
                inputMesh
            } );

            const relative = path.relative( inputDir, inputMesh );
            let destination = path.join( outputDir, relative );
            destination = destination.substring( 0, destination.length - 4 ) + "." + format;
            const destDir = path.dirname( destination );
            await desk.FileSystem.mkdirpAsync( destDir );

            await desk.Actions.executeAsync( {
                action : "copy",
                source : path.join( clean.outputDirectory, "mesh." + format ),
                destination
            } );

            cb();

        }, true );

    };

    ACVD.getColorFromValue = function( value, maxValue ) {

        value = Math.min( 1.0, value / maxValue );
        return chroma( 360 * ( 0.3 +  0.7 * value ), value , 1, "hsv" ).rgb();

    };

    ACVD.getDistanceBetweenMeshes = async function( sourceMesh, targetMesh ) {

        const distance =  await desk.Actions.executeAsync( {
            action : "M2MDistances", sourceMesh, targetMesh } );

        return distance.outputDirectory + "output.vtk";

    };

    ACVD.getHausdorffDistanceBetweenMeshes = async function( sourceMesh, targetMesh ) {

        const distance =  await desk.Actions.executeAsync( {
            action : "M2MHausdorff", sourceMesh, targetMesh, offsetSubdivisions : 1, stdout : true } );

        return distance.stdout.split( '\n' ).filter( l => l.startsWith( "Hausdorff" ) )
            .pop().split( " " ).map( parseFloat ).pop();

    };

    ACVD.getAverageDistanceBetweenMeshes = async function( sourceMesh, targetMesh ) {

        const distance =  await desk.Actions.executeAsync( {
            action : "M2MDistances", sourceMesh, targetMesh, stdout : true } );

        return parseFloat( distance.stdout.split( "\n" ).slice( -2 ).join( "" ).split( ":" ).pop().trim() );

    };

    ACVD.getColorBar = function( maxValue ) {

        const width = 1255, height = 1;
        const colorBar = new qx.ui.embed.Canvas().set( { syncDimension: false } );
        colorBar.set( { decorator :'tooltip', canvasHeight : height, canvasWidth : width } );
        const imageData = colorBar.getContext2d().createImageData( width, height );
        const data = imageData.data.fill( 255 );

        for ( let i = 0; i < width; i++ )
            data.set( this.getColorFromValue( i,  width ), 4 * i );

        colorBar.getContext2d().putImageData( imageData, 0, 0 );
        const minLabel = new qx.ui.basic.Label( "0" );
        const maxLabel = new qx.ui.basic.Label( "" + maxValue );
        const font = new qx.bom.Font( 50, [ "Arial" ] );
        for ( let label of [ minLabel, maxLabel ] ) label.setFont( font );

        const container = new qx.ui.container.Composite();
        container.setLayout( new qx.ui.layout.HBox( 5 ) );
        container.setUserData( "maxLabel", maxLabel );
        container.add( minLabel );
        container.add( colorBar, { flex : 1 } );
        container.add( maxLabel );
        return container;

    };

    ACVD.geometry2TextOBJ = function( geometry, convertToQuads = false ) {

        const lines = [];
        const pos = geometry.attributes.position;
        const indices = geometry.index.array;

        for (let i = 0; i < pos.count; i++ ) {

            const vertex = [ "v", pos.getX( i ), pos.getY( i ), pos.getZ( i ) ];
            lines.push( vertex );

        }

        if ( convertToQuads ) {

            for (let i = 0; i < indices.length / 6; i++ ) {

                const quad = [ "f" ];
                quad.push( indices[ 6 * i ] + 1 );
                quad.push( indices[ 6 * i + 1 ] + 1 );
                quad.push( indices[ 6 * i + 2 ] + 1 );
                quad.push( indices[ 6 * i + 5 ] + 1 );
                lines.push( quad );

            }

        } else {

            for (let i = 0; i < indices.length / 3; i++ ) {

                const tri = [ "f" ];
                tri.push( indices[ 3 * i ] + 1 );
                tri.push( indices[ 3 * i + 1 ] + 1 );
                tri.push( indices[ 3 * i + 2 ] + 1 );
                lines.push( tri );

            }

        }

        return lines.map( l => l.join( " " ) ).join( "\n" );

    };

    ACVD.downloadOBJ = function( mesh, fileName, convertToQuads = false ) {

        ACVD.downloadText( ACVD.geometry2TextOBJ( mesh.geometry, convertToQuads ), fileName );

    };

    ACVD.downloadText = function( text, fileName ) {

        const element = document.createElement('a');
        element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
        element.setAttribute('download', fileName);
        element.style.display = 'none';
        document.body.appendChild(element);
        element.click();
        document.body.removeChild(element);

    };


    ACVD.downloadArchive = async function( files, opts ) {

        const {

            tempDir = "data/ZipArchive",
            archiveName = "archive.zip",
            newNames = {}

        } = opts || {};

        if ( tempDir.split( "/" ).length < 2 ) throw( "cannot create temp dir from root");

        await desk.Actions.executeAsync( {

            action : "delete_directory",
            directory : tempDir

        } );

        await desk.Actions.executeAsync( {

            action : "mkdirp",
            directory : tempDir

        } );

        const join = window.require( "path" ).join;

        for ( let file of files ) {

            const outputFile = join( tempDir, newNames[ file ] || file.split( "/" ).pop() );

            await desk.Actions.executeAsync( {

                action : "copy",
                source : file,
                destination : outputFile

            } );

        }

        const outputZip = join( tempDir, archiveName );

        await desk.Actions.executeAsync( {

            action : "compress_to_zip",
            output_zip : outputZip,
            input_file_list : tempDir

        } );

        desk.FileSystem.downloadFile( outputZip );

    };


    ACVD.loadMeshWithColors = async function( file, maxDistance, opts ) {

        const mesh = await desk.THREE.Loader.loadMesh( file, opts );
        const nVertices = mesh.geometry.attributes.position.count;
        mesh.geometry.attributes.Distance.itemSize = 4;
        mesh.geometry.attributes.Distance.count = nVertices;
        const colors = new THREE.BufferAttribute( new Float32Array( nVertices * 3 ), 3 );

        for (let i = 0; i < nVertices; i++) {

            const distance = mesh.geometry.attributes.Distance.array[ 4 * i ];
            colors.set( this.getColorFromValue( distance, maxDistance ).map( v => v / 255 ), i * 3 );

        }

        mesh.geometry.setAttribute( 'color', colors );
        mesh.material.vertexColors = true;
        mesh.material.needsUpdate = true;
        return mesh;

    };

    ACVD.loadOBJ = async function ( file ) {

        const meshTxt = await desk.FileSystem.readFileAsync( file );
        const lines = meshTxt.split( "\n" ).map( l => l.split( " "  ) );
        const vertices = [];
        const indices = [];

        for ( let line of lines ) {
            switch ( line[ 0 ] ) {
                case "v" :
                    line.shift();
                    vertices.push( ...line.slice( 0, 3 ) );
                    break;

                case "f" :
                    line.shift();
                    line = line.map( n => n.split( "/" )[ 0 ] );
                    if ( line.length > 3 )
                        indices.push( line[ 0 ] - 1, line[ 1 ] - 1, line[ 2 ] - 1,
                            line[ 0 ] - 1, line[ 2 ] - 1, line[ 3 ] - 1);
                    else
                        indices.push( line[ 0 ] - 1, line[ 1 ] - 1, line[ 2 ] - 1 );
                    break;
                default:
            }
        }

        const geometry = new THREE.BufferGeometry();
        geometry.setIndex( indices );
        geometry.setAttribute( 'position', new THREE.BufferAttribute( new Float32Array( vertices ), 3 ) );
        geometry.computeVertexNormals();
        geometry.computeBoundingBox();
        return( geometry );

    };

    ACVD.readCSV = async function ( file ) {

        const txt = await desk.FileSystem.readFileAsync( file );
        return txt.split( ',' ).map( parseFloat ).filter( n => !isNaN( n ) );

    };


    window.ACVD = ACVD;
}
