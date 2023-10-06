"use strict";

{
    const { async, desk } = window;

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

    window.ACVD = ACVD;
}