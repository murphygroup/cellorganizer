Overview
------

1. Prerequisites for running CellOrganizer on Galaxy 
2. Downloading Pre-Trained Model on Galaxy
3. Synthesizing synthetic images and SBML Instance from Model
4. Uploading and viewing your realistic geometry on CellBlender
5. Creating a Lotka-Volterra Simulation with your geometry


Prerequisites
------
1. Galaxy ( 2 options)

   a. Locally installed version of `Galaxy <https://github.com/murphygroup/cellorganizer-galaxy-tools>`_ with Matlab R2018b
   b. Access to the public server `here <http://galaxy3.compbio.cs.cmu.edu:9000>`_

2. Installed version of `Blender with the CellBlender package <https://mcell.org/download.html>`_

Generating SBML Instance from Pretrained SPHARM-RPDM Model on Galaxy
------
1. Log in on Galaxy
2. From the right panel (History panel) click on the **“Gear Icon”**
3. From the drop down menu click on **“Create New”**

    .. figure:: galaxy-spharm_tutorial/images/create_new_history.png  

4. To rename the history. Double click on **“Unnamed history”** and rename it to **“SPHARM Model”**. Then click enter.
    
    .. figure:: galaxy-spharm_tutorial/images/renaming_history1.png


5. Annotate history. Click on the tag icon to add tags to this history. Add **“train”** and **“vesicle”** as tags. Click enter after each tag.

   .. figure:: galaxy-spharm_tutorial/images/add_tags.png

6. Import dataset from shared data. Click on **“Shared Data”** at the top of the screen then click on **“Data Libraries”** from the dropdown menu. 
    
    .. figure:: galaxy-spharm_tutorial/images/data_libraries.png

7. In the following page, follow the links for:

   a. Click **“Data Libraries”** from the drop down menu.
   b. Click on **“Generative Models”**
   c. Click on **“HeLa”**
   d.  Click on **“3D”**
   e.  Click on **“3D_HeLa_LAMP2_SPHARM_vesicle_model.mat”**. 

8.  Then ticking off the box to the left of its name, then click on **“To History”** and select the history called **“SPHARM Model”**

    .. figure:: galaxy-spharm_tutorial/images/to_history.png

9. Then click on **“to history”** button in the top menu. A dialog window with appear with the current history on board pre-selected or you can create a new one as it gives you that option as well.
   * Click on the **“import”** button if your history is already pre-selected this will import both datasets into your history. Once the images are imported **a green box in the top right corner** will appear, click on it so it will take you to the history with the images imported

10. Or you can also click on the |galaxy_button| icon in the top left corner of the screen  to return to the home page. 

    .. |galaxy_button| image:: galaxy-spharm_tutorial/images/galaxy_button.png  

11. Then, click on **“Synthesis”** under the **“Tools”** menu, and follow the link to **“Synthesize an instance from multiple models trained in CellOrganizer”**
    
    .. figure:: galaxy-spharm_tutorial/images/Tools_panel.png

12. Select the model from the list and select **"Synthesize from all models"** as the synthesis option.

    .. figure:: galaxy-spharm_tutorial/images/Synthesis_from_all_models.png

13. To save the output as an image and SBML mesh instance, click the YES button under Output Options for: OMETIFF, SBML Spatial 3, and Indexed Image

    .. figure:: galaxy-spharm_tutorial/images/Ome_tiff_options.png

14. In the **“Advanced Options”**, match the following image:
    
    .. figure:: galaxy-spharm_tutorial/images/adv_options.png

15. Once all the information is complete click on the |execute| button, that will close the options panel. A green box will be displayed indicating that the demo is being run and a new item in the history will be added with the model ran. 

    .. |execute| image:: galaxy-spharm_tutorial/images/execute_button.png
    
    * You should see your generated outputs in the right sidebar
        .. figure:: galaxy-spharm_tutorial/images/outputs1_right_sidebar.png
    
16. You can view the indexed image by clicking the eye icon next to the name
    
    .. figure:: galaxy-spharm_tutorial/images/view_result_right_sidebar.png

Importing Generated SBML instance into CellBlender
------

1. Download the SBML instance from Galaxy clicking the eye icon

    .. figure:: galaxy-spharm_tutorial/images/SBML_Galaxy.png

2. Next, open up Blender with CellBlender pre-installed. Initialize CellBlender.
     
    .. figure:: galaxy-spharm_tutorial/images/initialize_blender.png

3. Import the downloaded SBML instance by going to: **File > Import > BioNetGen/SBML Model(.bng, ./xml)**.  You should now see your imported SBML instance. Use the scroll-pad and mouse to move around and investigate the geometry.

    .. figure:: galaxy-spharm_tutorial/images/Import_blender.png

Create a Lotka-Volterra Simulation with our realistic geometry
------
1. Next step is to then import a .txt file, located at XXXXX, that includes the preset reactions for our simulation. Go to: **File >Import >CellBlender Model(text/pickle)**

    .. figure:: galaxy-spharm_tutorial/images/SBML_instance.png

2. Next, we have to rescale and color our simulated particles. Under the **"Molecules"** button, set the scale of both **"prey"** and **"predator"** to 20.0. Set the color of **"prey"** to blue and **"predator"** to red. 
    
    .. image:: galaxy-spharm_tutorial/images/color_properties_CB1.png
        :width: 49 %

    .. image:: galaxy-spharm_tutorial/images/color_properties_CB2.png
        :width: 49 %
3. Then, save the file as SPHARM_Model_Sim.blend. Next, you should see the Run button appear under the Run Simulation tab. Click that.

    .. figure:: galaxy-spharm_tutorial/images/run_simulation_CB.png
    
    **Note:** It's possible that the Run button doesn't appear. This may be caused by the Mcell binary path not being set if not by default. Go to the Preferences tab under CellBlender and navigate to the option to set Mcell Binary Path. Depending on your device, this path should then be set to:
    
    * Linux: `/home/[user]/blender-[version]/[version number]/scripts/addons/cellblender/bin/mcell`  
    * Windows: `C:\Users\[user]\AppData\Romaing\Blender Foundation\Blender\[user]\scripts\addons\cellblender\bin\mcell`  
    * `C:\ProgramData\Blender Foundation\Blender\[user folder]\scripts\addons\cellblender\bin\mcell`
    
    with [user] and [version number] depending on your device.

4. This should produce a simulation similar to the one shown:
    
    .. figure:: galaxy-spharm_tutorial/images/CellBlender_FullScreen_gif.gif


