
const defaultParameters = [
	{
		name: "inputMesh",
		alias: "input_mesh",
		type: "file",
		required: true
	},
	{
		name: "numberOfVertices",
		alias: "number_of_desired_vertices",
		info: "number of desired vertices for the output mesh",
		type: "int",
		min: 0,
		required: true
	},
	{
		name: "gradation",
		type: "float",
		info: "set to 0 for uniform meshing, increase (example : 2.5) to add curvature-adapted meshing",
		required: true,
		min: 0,
		max: 30,
		defaultValue: 0
	},
	{
		prefix: "-m ",
		info: "force the output mesh to be manifold",
		name: "forceManifold",
		alias: "force_manifold",
		type: "int",
		min: 0,
		max: 1
	},
	{
		prefix : "-s ",
		info : "subsampling threshold : subdivide the mesh until its number of points is large enough i.e. numberOfInputVertices / number_of_desired_vertices > subsamplingThreshold",
		name : "subsamplingThreshold",
		type : "int"
	},
	{
		prefix : "-l ",
		info : "split input edges longer than ( averageLength * lengthRatio )",
		name : "lengthRatio",
		type : "float"
	},
	{
		text: "-d 0"
	}
];

const QParameters = [
	{
		prefix: "-q ",
		info: "quadrics level (0 : no quadrics, 3 : full quadrics)",
		name: "quadricsLevel",
		type: "int",
		min: 0,
		max: 3,
		defaultValue: 3
	}
]

const opts = {
    actions : {
		acvd: {
			description: "Performs simplification/remeshing of 3D triangular meshes",
			parameters : defaultParameters,
			executable: "bin/ACVD"
		},
		acvdp: {
			description: "Performs simplification/remeshing of 3D triangular meshes (parallel version)",
			parameters: defaultParameters,
			executable: "bin/ACVDP"
		},
		acvdq: {
			description: "Performs simplification/remeshing of 3D triangular meshes (usin QEM, more accurate but slower)",
			parameters: [ ...defaultParameters, ...QParameters ],
			executable: "bin/ACVDQ"
		},
		acvdqp: {
			description: "Performs simplification/remeshing of 3D triangular meshes (usin QEM, more accurate but slower, parallel version)",
			parameters: [ ...defaultParameters, ...QParameters ],
			executable: "bin/ACVDQP"
		},
    }
}

module.exports = opts;
