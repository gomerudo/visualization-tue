/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package gui;

import javax.swing.JOptionPane;
import util.ControlOptions;
import volvis.RaycastRenderer;

/**
 *
 * @author michel
 */
public class RaycastRendererPanel extends javax.swing.JPanel {

    RaycastRenderer renderer;
    TransferFunctionEditor tfEditor = null;
    TransferFunction2DEditor tfEditor2D = null;
    
    /**
     * Creates new form RaycastRendererPanel
     */
    public RaycastRendererPanel(RaycastRenderer renderer) {
        initComponents();
        this.renderer = renderer;
    }

    public void setSpeedLabel(String text) {
        renderingSpeedLabel.setText(text);
    }
    
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        buttonGroup1 = new javax.swing.ButtonGroup();
        jLabel1 = new javax.swing.JLabel();
        renderingSpeedLabel = new javax.swing.JLabel();
        slicerButton = new javax.swing.JRadioButton();
        mipButton = new javax.swing.JRadioButton();
        compositingButton = new javax.swing.JRadioButton();
        tf2dButton = new javax.swing.JRadioButton();
        shadingCheckbox = new javax.swing.JCheckBox();
        nSlices = new javax.swing.JTextField();
        jLabel2 = new javax.swing.JLabel();
        lowResolution = new javax.swing.JCheckBox();

        jLabel1.setText("Rendering time (ms):");

        renderingSpeedLabel.setText("jLabel2");

        buttonGroup1.add(slicerButton);
        slicerButton.setSelected(true);
        slicerButton.setText("Slicer");
        slicerButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                slicerButtonActionPerformed(evt);
            }
        });

        buttonGroup1.add(mipButton);
        mipButton.setText("MIP");
        mipButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                mipButtonActionPerformed(evt);
            }
        });

        buttonGroup1.add(compositingButton);
        compositingButton.setText("Compositing");
        compositingButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                compositingButtonActionPerformed(evt);
            }
        });

        buttonGroup1.add(tf2dButton);
        tf2dButton.setText("2D Transfer function");
        tf2dButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                tf2dButtonActionPerformed(evt);
            }
        });

        shadingCheckbox.setText("Volume shading");
        shadingCheckbox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                shadingCheckboxActionPerformed(evt);
            }
        });

        nSlices.setText(ControlOptions.N_SLICES + "");
        nSlices.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                nSlicesActionPerformed(evt);
            }
        });

        jLabel2.setText("Slices:");

        lowResolution.setText("Low resolution");
        lowResolution.setVisible(false);
        lowResolution.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                lowResolutionActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel1)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(renderingSpeedLabel))
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addComponent(lowResolution)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jLabel2)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(nSlices))
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(compositingButton)
                                .addComponent(tf2dButton)
                                .addComponent(mipButton)
                                .addComponent(slicerButton)
                                .addComponent(shadingCheckbox)))))
                .addContainerGap(339, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(renderingSpeedLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(slicerButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(mipButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(compositingButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(tf2dButton)
                .addGap(18, 18, 18)
                .addComponent(shadingCheckbox)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(nSlices, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel2))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(lowResolution)
                .addContainerGap(122, Short.MAX_VALUE))
        );
    }// </editor-fold>//GEN-END:initComponents

    private void mipButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_mipButtonActionPerformed
        ControlOptions.VISUALIZATION_MODE = ControlOptions.MIP_OPT;
        renderer.changed();
    }//GEN-LAST:event_mipButtonActionPerformed

    private void slicerButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_slicerButtonActionPerformed
        ControlOptions.VISUALIZATION_MODE = ControlOptions.SLICER_OPT;
        renderer.changed();
    }//GEN-LAST:event_slicerButtonActionPerformed

    private void compositingButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_compositingButtonActionPerformed
        ControlOptions.VISUALIZATION_MODE = ControlOptions.COMPOSITING_OPT;
        renderer.changed();
    }//GEN-LAST:event_compositingButtonActionPerformed

    private void tf2dButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_tf2dButtonActionPerformed
        ControlOptions.VISUALIZATION_MODE = ControlOptions.TF2D_OPT;
        renderer.changed();
    }//GEN-LAST:event_tf2dButtonActionPerformed

    private void shadingCheckboxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_shadingCheckboxActionPerformed
        
        ControlOptions.SHADDING = shadingCheckbox.isSelected();
        renderer.changed();
    }//GEN-LAST:event_shadingCheckboxActionPerformed

    private void nSlicesActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_nSlicesActionPerformed
        try {
            int value = Integer.parseInt(nSlices.getText());
            
            ControlOptions.N_SLICES = value;
        } catch (NumberFormatException e) {
        }
        renderer.changed();
    }//GEN-LAST:event_nSlicesActionPerformed

    private void lowResolutionActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_lowResolutionActionPerformed
        ControlOptions.LOW_RESOLUTION = lowResolution.isSelected();
        renderer.changed();
    }//GEN-LAST:event_lowResolutionActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.ButtonGroup buttonGroup1;
    private javax.swing.JRadioButton compositingButton;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JCheckBox lowResolution;
    private javax.swing.JRadioButton mipButton;
    private javax.swing.JTextField nSlices;
    private javax.swing.JLabel renderingSpeedLabel;
    private javax.swing.JCheckBox shadingCheckbox;
    private javax.swing.JRadioButton slicerButton;
    private javax.swing.JRadioButton tf2dButton;
    // End of variables declaration//GEN-END:variables
}
