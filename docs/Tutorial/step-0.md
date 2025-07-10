## Write a simple input file:

In your input file (`*.json` file), the layout of the block should resemble the following:
```
"mesh":{
		"type":"asfem",
		"dim":1,
		"nx":20,
		"xmax":1.5,
		"meshtype":"edge2",
		"savemesh":true
	}
```
In this context, `"mesh"` represents the function or module we want to use in AsFem, and it defines a single block in the JSON file. Within this block, you can include multiple options for different purposes.

In the aforementioned example, `"mesh"` defines a block for mesh generation, and within the `{}` bracket pair, you can specify the mesh type (`"type":"asfem"`), the dimension of the mesh (`"dim":1`), and other relevant details.


To create a minimal input file for AsFem, you will need to include at least the following blocks:
```
"mesh": to define the finite element mesh
"dofs": to define the degree of freedom
"elements": to define the PDE/ODE
"job": to define the static/transient analysis job
```


Within the `"dofs"` block, you can specify the name or list of names of the *Degrees of Freedom* (**DoFs**) that will be displayed in Paraview when you open the result file (`vtu` file). Furthermore, when applying boundary conditions, you will need to provide the name of the relevant **DoFs**.
```
"dofs":{
		"names":["phi"]
	}
```

and the element or module you'd like to use:

Additionally, you will need to specify the element or module you would like to use in the `"elements"` block. This block defines the residual and jacobian for your PDE/ODE.
```
"elements":{
		"elmt1":{
			"type":"laplace",
			"dofs":["phi"],
			"material":{
				"type":"constpoisson",
				"parameters":{
					"sigma":1.0,
					"f":1.0e1
				}
			}
		}
	}
```


In the end, you will also need to provide the necessary information for the analysis within the `"job"` block.
```
"job":{
		"type":"static",
		"print":"dep",
		"restart":true
	}
```

To summarize, in order to create a basic input file for AsFem, you will need to include the following blocks:
```
{
	"mesh":{
		"type":"asfem",
		"dim":2,
		"nx":50,
		"ny":50,
		"meshtype":"quad4",
		"savemesh":true
	},
	"dofs":{
		"names":["phi"]
	},
	"elements":{
		"elmt1":{
			"type":"laplace",
			"dofs":["phi"],
			"material":{
				"type":"constpoisson",
				"parameters":{
					"sigma":1.0,
					"f":1.0e1
				}
			}
		}
	},
	"job":{
		"type":"static",
		"print":"dep"
	}
}
```
In addition to performing a full FEM analysis workflow, you can use the `--read-only` option with AsFem to parse and validate the input file without executing the simulation. This option allows you to check the syntax and structure of the input file without actually running a full analysis.

You can find examples of all the possible blocks that can be included in an AsFem input file in the `test_input/input-template.json` file.

In the following steps, we will provide detailed explanations of each block and their functions in AsFem. You can find example input files in the `examples/tutorial` folder to help you understand the usage of each block.
