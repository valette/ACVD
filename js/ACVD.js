"use strict";

{
    const { require, async, desk } = window;

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

    window.ACVD = ACVD;
}
